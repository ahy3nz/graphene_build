import os
import numpy as np
import mdtraj
import itertools
from scipy.spatial import distance
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import subprocess

constants = ['50', '500']
trials = ['a', 'b', 'c']
curr_dir = os.getcwd()
for combo in itertools.product(constants,  trials):
    thing = 'k{}_{}'.format(combo[0], combo[1])
    os.chdir(os.path.join(curr_dir, thing))
    print(thing)
    if not os.path.isfile('pull_nopbc.xtc'):
        p = subprocess.Popen('echo 0 | gmx trjconv -f pull.xtc -s pull.tpr -pbc mol -o pull_nopbc.xtc', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()
    traj = mdtraj.load('pull_nopbc.xtc', top='em.gro')
    box_dims = traj.unitcell_lengths[0]
    pullf = np.loadtxt('pull_pullf.xvg', comments=['@', '#'])
    raw_force = pullf
    
    
    # Reduce amount of data in pullf to match len of traj
    block_size = 10
    remainder = pullf.shape[0] % block_size
    pullf = pullf[0:-remainder,1]
    n_blocks = int(pullf.shape[0]/block_size)
    pullf = pullf.reshape(n_blocks, block_size)
    pullf = np.mean(pullf, axis=1)
    
    graphene_atom = int(open('index.ndx', 'r').readlines()[-1]) - 1
    integrand = []
    for i, force in enumerate(pullf):
        dx = distance.euclidean(traj.xyz[i, graphene_atom, :], 
                                traj.xyz[i+1, graphene_atom, :])
        # Note the square box
        if dx > box_dims[0]:
            dx -= box_dims[0]

        integrand.append(force*dx)
    
    work_profile = np.cumsum(integrand)
    all_times = traj.time[1:]
    
    fig, ax = plt.subplots(1,1)
    ax.plot(raw_force[:,0], raw_force[:,1], color='r')
    ax.set_ylabel(r"Force ($\frac{kJ}{mol*nm}$)", color='r')
    for tick in ax.get_yticklabels():
        tick.set_color('r')
    
    ax2 = ax.twinx()
    ax2.plot(all_times, work_profile, color='b')
    np.savetxt('work_profile.dat', np.column_stack((all_times, work_profile)))
    ax2.set_ylabel(r"Work ($\frac{kJ}{mol}$)", color='b')
    for tick in ax2.get_yticklabels():
        tick.set_color('b')
    
    ax.set_xlabel("Time (ps)")
    fig.tight_layout()
    fig.savefig("profile.png", transparent=True)
    plt.close(fig)
    
