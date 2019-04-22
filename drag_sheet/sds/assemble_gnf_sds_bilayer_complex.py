import os
import shutil
import glob
import subprocess

import load_graphene
#sds_vals = [10, 20, 30, 40, 50, 60] # number of SDS on both sides of graphene
sds_vals=[70,80]
angles = ['0']
trials = ['a', 'b', 'c']

PATH_TO_MOLECULES = '/raid6/homes/ahy3nz/Programs/graphene_build/single_molecules/'

PATH_TO_EM='/raid6/homes/ahy3nz/Trajectories/charmm/sds/inputs/'
curr_dir = os.getcwd()
for n_sds in sds_vals:
    os.chdir(curr_dir)
    new_folder = '{}sds'.format(str(n_sds))
    if not os.path.isdir(new_folder):
        os.makedirs(new_folder)
    os.chdir(os.path.join(curr_dir, new_folder))
    for trial in trials:
        os.chdir(os.path.join(curr_dir, new_folder))
        out_folder = "k50_0_{}".format(trial)
        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)
        os.chdir(os.path.join(curr_dir, new_folder, out_folder))
        # Create bilayer, gnf, sds + ion system
        # Write pulling mdp
        print(os.getcwd())
        print("constructing ...")
        load_graphene.construct_system(n_load_per_side=int(n_sds/2),
                dist_from_sheet=0.2)

        # add water
        print("solvating ...")
        solvate_command = 'gmx solvate -cp composite_rotated.gro -cs {}/water.gro -p topol.top -o solvated_{}sds.gro -scale 0.57'.format(PATH_TO_MOLECULES,n_sds)
        p = subprocess.Popen(solvate_command, shell=True, 
                stdout=open('solv.log', 'w'), 
                stderr=open('solv.err', 'w'))
        p.wait()

        # create index
        print("indexing ...")
        index_command = 'echo q | gmx make_ndx -f solvated_{}sds.gro'.format(n_sds)
        p = subprocess.Popen(index_command, shell=True,
                stdout=open('index.log', 'w'),
                stderr=open('index.err', 'w'))
        p.wait()
        with open('index.ndx', 'a') as f:
            f.write("\n[ anchor2] \n")
            f.write('1007\n')

        # run EM
        print("em-ing ...")
        em_command = 'gmx grompp -f {}/em.mdp -c solvated_{}sds.gro -p topol.top -o em'.format(PATH_TO_EM, n_sds)
        p = subprocess.Popen(em_command, shell=True, 
                stdout=open('em_grompp.log', 'w'),
                stderr=open('em_grompp.err', 'w'))
        p.wait()
        em_command = 'gmx mdrun -deffnm em'
        p = subprocess.Popen(em_command, shell=True, 
                stdout=open('em_run.log', 'w'),
                stderr=open('em_run.err', 'w'))
        p.wait()

        # grompp for angled insertion
        print("grompping ...")
        pull_command = 'gmx grompp -f angled_insertion.mdp -c em.gro -p topol.top -n index.ndx -maxwarn 1 -o pull'
        p = subprocess.Popen(pull_command, shell=True, 
                stdout=open('pull_grompp.log', 'w'),
                stderr=open('pull_grompp.err', 'w'))
        p.wait()

        # Reorganize and move files
        #to_delete = glob.glob('#*')
        #for foo in to_delete:
        #    os.unlink(foo)
        #
        #extensions = ['*.gro', '*.log', '*.edr', '*.top', '*.ndx', 
        #    '*.tpr', '*.xvg', '*.trr', '*.mdp']
        #to_move = []
        #for ext in extensions:
        #    to_move.extend(glob.glob(ext))
        #for foo in to_move:
        #    shutil.move(foo, os.path.join(out_folder, foo))
