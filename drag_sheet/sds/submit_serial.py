import os
import itertools as it

import script_utils
from operations import submit_job, compute_work, compute_penetration
from multiprocessing import Pool

def hello(arg0):
    print((arg0))

n_nodes = 2
curr_dir = os.getcwd()
sds_folders = ['10sds', '20sds', '30sds', '40sds', '50sds', '60sds']
k_folders = ['k50']
angle_folders = ['0']
trials = ['a', 'b', 'c']
all_trials = [*it.product(sds_folders, k_folders, angle_folders,trials)]
all_sims = [os.path.join(curr_dir, '{0}/{1}_{2}_{3}'.format(*combo)) for combo in all_trials]
#print(all_sims)
iterator = iter(all_sims)
with Pool(n_nodes) as p:
    #p.map(compute_work, all_sims)
    p.map(compute_penetration, all_sims)
#for thing in zip(iterator, iterator):
#    print(thing)
#    import pdb; pdb.set_trace()
#    with Pool(n_nodes) as p:
#        p.starmap(hello, thing)
    
#for combo in it.product(sds_folders, k_folders, angle_folders,trials):
#    sim_dir = os.path.join(curr_dir, '{0}/{1}_{2}_{3}'.format(*combo))
#    print("{}".format(sim_dir))
#    compute_work(sim_dir)
#    compute_penetration(sim_dir)

