import pandas as pd
import numpy as np
import os
import itertools as it
import json

df = pd.DataFrame()
curr_dir = os.getcwd()
sds_folders = ['10sds', '20sds', '30sds', '40sds', '50sds', '60sds']
k_folders = ['k50']
angle_folders = ['0']
trials = ['a', 'b', 'c']
all_statepoints = [*it.product(sds_folders, k_folders, angle_folders)]
#all_sims = [os.path.join(curr_dir, '{0}/{1}_{2}_{3}'.format(*combo)) for combo in all_trials]

for statepoint in all_statepoints:
    coating = statepoint[0]
    k = statepoint[1]
    angle = statepoint[2]
    all_pens = []
    all_pen_stds = []
    all_works = []
    all_work_stds = []
    for trial in trials:
        chdir = os.path.join(curr_dir, 
                '{0}/{1}_{2}_{3}'.format(*statepoint, trial))
        
        os.chdir(chdir)
        print(chdir)
        work_profile = np.loadtxt('work_distance_profile.dat')
        all_works.append(work_profile[-1, 1])
        all_work_stds.append(np.std(work_profile[-100:, 1]))
        penetration = json.load(open('penetration.json', 'r'))
        all_pens.append(penetration['distance'])
        all_pen_stds.append(penetration['err'])
    to_add = {'work': np.mean(all_works),
            'work_std': np.sqrt(np.sum([a**2 for a in all_work_stds])),
            'penetration': np.mean(all_pens),
            'penetration_std': np.sqrt(np.sum([a**2 for a in all_pen_stds])),
            'k':k,
            'angle': angle,
            'coating': coating
            }
    df = df.append(to_add, ignore_index=True)
            
columns=['coating', 'k', 'angle', 'work', 'work_std',
                'penetration', 'penetration_std']
df = df[columns]
os.chdir(curr_dir)
df.to_csv('results.csv')
