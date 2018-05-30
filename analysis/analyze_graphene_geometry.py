import numpy as np
import mbuild as mb
import mdtraj

############################
## Compute the angle of atack between graphene sheet and bilayer
## Note the hardcoded indices for frame and atom
######################

#cmpd = mb.load('run.gro')
traj = mdtraj.load('run.gro')
point1 = traj.xyz[0,49568,:]
point2 = traj.xyz[0,50575,:]
cmpd_vector = point2 - point1
ref_vector = [0,0,-1]


theta = mb.coordinate_transform.angle(cmpd_vector, ref_vector)
print(np.rad2deg(theta))

