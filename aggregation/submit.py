import os
import operations
import script_utils

folders = ['trial4', 'trial5', 'trial6']
main_dir = os.getcwd()
jobids = [None] * 2
for i, folder in enumerate(folders):
    os.chdir(os.path.join(main_dir, folder))
    body = 'cd {}\n'.format(os.getcwd())
    body += "gmx mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm npt"
    with open('submit.pbs', 'w') as f:
        script_utils.write_rahman_script(f,jobname="graphene_agg_{}".format(folder),
                body=body) 
    operations.submit_job('submit.pbs', jobids, 2, i)
