import os
import glob
import subprocess
import shutil

import mbuild as mb

import scripts.bilayer as bilayer

##################
## Using gromacs utilities to build the system
## But using mbuild functions to write the topology
##################

header = """; Include forcefield parameters
#include "/raid6/homes/ahy3nz/Trajectories/graphene_studies/graphene_files/force-field/ffmix.itp"
#include "/raid6/homes/ahy3nz/Trajectories/graphene_studies/graphene_files/force-field/molecules.itp"
#include "/raid6/homes/ahy3nz/Trajectories/graphene_studies/graphene_files/force-field/ions.itp" """

ingredients_folder = 'ingredients'
out_folder = 'trial1'

box = [10, 10, 10]
cmpd = mb.Compound()
sheet = 'sheet.gro'
components = [ ('GRA_5', 1, 'sheet.gro', mb.load('{}/sheet.gro'.format(ingredients_folder))), 
                ('DOPC', 100, 'dopc.gro', mb.load('{}/dopc.gro'.format(ingredients_folder))),
                ('CHL1', 100, 'chl1.gro', mb.load('{}/chl1.gro'.format(ingredients_folder))),
                ('SOL', 1000, 'SOL.gro', mb.load('{}/SOL.gro'.format(ingredients_folder))) ]
outfile = '{}/bulk.gro'.format(out_folder)
outtop = '{}/bulk.top'.format(out_folder)
old_files = glob.glob('*bulk.gro*')
for thing in old_files:
    os.remove(thing)
with open('{}/build.out'.format(out_folder),'w') as stdoutfile, open('{}/build.err'.format(out_folder), 'w') as stderrfile:
    for i, (name, nmol, grofile, cmpd_i)  in enumerate(components):
        if i == 0:
            cmd = ('gmx insert-molecules -ci {4}/sheet.gro -nmol 1 -o {0} -box '
                '{1} {2} {3}'.format(outfile, *box, ingredients_folder))
        else:
            cmd = ('gmx insert-molecules -f {0} -ci {6}/{1} '
                '-nmol {2} -o {0} -box {3} {4} {5}'.format(outfile, grofile, 
                    nmol, *box, ingredients_folder))
        print(cmd)
        p = subprocess.Popen(cmd,
                stdout=stdoutfile, stderr=stderrfile, shell=True)
        p.wait()
        for _ in range(nmol):
            to_add = mb.clone(cmpd_i)
            to_add.name = name
            cmpd.add(to_add)

    bilayer.write_gmx_topology(cmpd, outtop, header=header)


with open('{0}/em.out'.format(out_folder), 'w') as log:
    p = subprocess.Popen('gmx grompp -f mdp/em.mdp -c {0} '
            '-p {1} -o {2}/em'.format(outfile, outtop, out_folder),
            shell=True, stdout=log, stderr=log)
    p.wait()
    p = subprocess.Popen('gmx mdrun -deffnm {0}/em'.format(out_folder), 
            shell=True, stdout=log, stderr=log)
    p.wait()

with open('{0}/npt.out'.format(out_folder), 'w') as log:
    p = subprocess.Popen('gmx grompp -f mdp/npt.mdp -c {1}/em.gro '
            '-p {0} -o {1}/npt -maxwarn 1'.format(outtop, out_folder),
            shell=True, stdout=log, stderr=log)
    p.wait()
