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


box = [10, 10, 10]
cmpd = mb.Compound()
sheet = 'sheet.gro'
components = [ ('GRA_5', 1, 'sheet.gro', mb.load('sheet.gro')), 
                ('DOPC', 100, 'dopc.gro', mb.load('dopc.gro')),
                ('CHL1', 100, 'chl1.gro', mb.load('chl1.gro')),
                ('SOL', 1000, 'SOL.gro', mb.load('SOL.gro')) ]
outfile = 'bulk.gro'
old_files = glob.glob('*bulk.gro*')
for thing in old_files:
    os.remove(thing)
with open('std.out','w') as stdoutfile, open('std.err', 'w') as stderrfile:
    for i, (name, nmol, grofile, cmpd_i)  in enumerate(components):
        if i == 0:
            cmd = ('gmx insert-molecules -ci sheet.gro -nmol 1 -o {0} -box '
                '{1} {2} {3}'.format(outfile, *box))
        else:
            cmd = ('gmx insert-molecules -f {0} -ci {1} '
                '-nmol {2} -o {0} -box {3} {4} {5}'.format(outfile, grofile, 
                    nmol, *box))
        print(cmd)
        p = subprocess.Popen(cmd,
                stdout=stdoutfile, stderr=stderrfile, shell=True)
        p.wait()
        for _ in range(nmol):
            to_add = mb.clone(cmpd_i)
            to_add.name = name
            cmpd.add(to_add)

bilayer.write_gmx_topology(cmpd, 'bulk.top', header=header)
