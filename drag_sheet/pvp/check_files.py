import os
import glob
import itertools as it

curr_dir = os.getcwd()
sds_folders = ['10pvp', '20pvp', '30pvp']
k_folders = ['k50']
angle_folders = ['0']
trials = ['a', 'b', 'c']
for combo in it.product(sds_folders, k_folders, angle_folders,trials):
    sim_dir = os.path.join(curr_dir, '{0}/{1}_{2}_{3}'.format(*combo))
    if os.path.isdir(sim_dir):
        os.chdir(sim_dir)
        to_remove = glob.glob('pull*')
        for thing in to_remove:
            if 'pull.tpr' not in thing:
                os.unlink(thing)
        #if not os.path.isfile('pull.tpr'):
            #print('{} error'.format(sim_dir))
    else:
        print("{} doesn't exist".format(sim_dir))
