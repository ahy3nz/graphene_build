import os
import itertools as it

curr_dir = os.getcwd()
sds_folders = ['10sds', '20sds', '30sds', '40sds', '50sds', '60sds', 
            '70sds', '80sds', '90sds', '100sds']
k_folders = ['k50']
angle_folders = ['0']
trials = ['a', 'b', 'c']
for combo in it.product(sds_folders, k_folders, angle_folders,trials):
    sim_dir = os.path.join(curr_dir, '{0}/{1}_{2}_{3}'.format(*combo))
    if os.path.isdir(sim_dir):
        os.chdir(sim_dir)
        if not os.path.isfile('pull.gro'):
            print(sim_dir)
    else:
        print("{} doesn't exist".format(sim_dir))
