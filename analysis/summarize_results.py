import os
import numpy as np 
import pandas as pd 
import pdb
trials = ['a', 'b','c']
constants = ['50', '100', '125']
#constants = ['250', '500', '1000']
angles = ['0']
#angles = ['0', '15', '30', '45']
df = pd.DataFrame()
curr_dir = os.getcwd()
for k in constants:
    for a in angles:
        all_forces = []
        all_works = []
        for trial in trials:
            name = 'k{}_{}_{}'.format(k, a, trial)
            os.chdir(os.path.join(curr_dir, name))
            forces = np.loadtxt('pull_pullf.xvg', comments=['@', '#'])
            max_force = np.max(forces[:,1])
            work = np.loadtxt('work_profile.dat')
            max_work = np.mean(work[int(0.9*work.shape[0]) : , 1])
            all_forces.append(max_force)
            all_works.append(max_work)
        to_add = {'k': [k], 'angle': [a], 
                'max_force': [np.mean(all_forces)], 'max_force_std': [np.std(all_forces)],
                'max_work': [np.mean(all_works)], 'max_work_std': [np.std(all_works)]}
        df = df.append(pd.DataFrame.from_dict(to_add))
os.chdir(curr_dir)
df.to_csv('summary_weak.csv')

