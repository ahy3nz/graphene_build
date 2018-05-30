import numpy as np
import scipy.integrate
import pdb

#################
## Basic script to look at forces and positions from gromacs simulation
## and integrate to get the work 
#############

pullf = np.loadtxt('pull_pullf.xvg', comments=['@', '#'])
pullx = np.loadtxt('pull_pullx.xvg', comments=['@', '#'])

first_integral = scipy.integrate.cumtrapz(pullf[:,1], pullx[:,1])
second_integral = scipy.integrate.cumtrapz(pullf[:,2], pullx[:,2])
total_work = first_integral + second_integral
total_work_new = [val for val in total_work]
total_work_new.append(total_work_new[-1])
total_force = -1*(pullf[:,1] + pullf[:,2])
print(total_force[-1], total_work_new[-1])
np.savetxt('pull_physics.dat', np.column_stack((pullx[:,0], total_force, total_work_new)))
