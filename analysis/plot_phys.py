import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

###################
## Script to plot force and work on same plot with twin axes
##############

data = np.loadtxt('pull_physics.dat')
fig, ax = plt.subplots(1,1, figsize=(8,5))
ax.plot(data[:,0], data[:,1], label='Force', color='red')
ax2 = ax.twinx()
ax2.plot(data[:,0], data[:,2], label='Work', color='blue')
ax.set_ylabel("Force", color='red', fontsize=20)
ax2.set_ylabel("Work", color='blue', fontsize=20)
ax.set_xlabel("Time (ps)")
fig.savefig('physics.jpg')
