import os
import numpy as np 
import pandas as pd 
import pdb
trials = ['a', 'b','c']
constants = ['50', '500']

df = pd.DataFrame()
curr_dir = os.getcwd()
for k in constants:
    all_forces = []
    all_works = []
    for trial in trials:
        name = 'k{}_{}'.format(k, trial)
        os.chdir(os.path.join(curr_dir, name))
        forces = np.loadtxt('pull_pullf.xvg', comments=['@', '#'])
        max_force = np.max(forces[:,1])
        work = np.loadtxt('work_profile.dat')
        #max_work = np.mean(work[int(0.9*work.shape[0]) : , 1])
        max_work = work[-1,1] 
        all_forces.append(max_force)
        all_works.append(max_work)
    to_add = {'k': [k], 
            'max_force': [np.mean(all_forces)], 'max_force_std': [np.std(all_forces)],
            'max_work': [np.mean(all_works)], 'max_work_std': [np.std(all_works)],
            'dna':'2dna'}
    df = df.append(pd.DataFrame.from_dict(to_add))
os.chdir(curr_dir)
df.to_csv('work.csv')

