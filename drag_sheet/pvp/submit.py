import os
import itertools as it

import script_utils
from operations import submit_job, compute_work, compute_penetration
from multiprocessing import Pool

n_nodes = 2
jobid = [None]*n_nodes
i=0
curr_dir = os.getcwd()
sds_folders = ['40pvp']
k_folders = ['k50']
angle_folders = ['0']
trials = ['a', 'b', 'c']
for combo in it.product(sds_folders, k_folders, angle_folders,trials):
    sim_dir = os.path.join(curr_dir, '{0}/{1}_{2}_{3}'.format(*combo))
    print("{}".format(sim_dir))
    if os.path.isdir(sim_dir):
        os.chdir(sim_dir)
        if os.path.isfile('pull.tpr') and not os.path.isfile('pull.gro'):
            cmd = 'cd {}\n'.format(sim_dir)
            cmd += 'gmx mdrun -deffnm pull -ntomp 8 -ntmpi 2 -gpu_id 01 -append -cpi pull_prev.cpt -pf pull_pullf.xvg -px pull_pullx.xvg\n'
            with open('pull.pbs', 'w') as f:
                script_utils.write_rahman_script(f, jobname='{0}_{3}'.format(*combo),
                        body=cmd)
            submit_job('pull.pbs',jobid,n_nodes,i)
            i += 1
