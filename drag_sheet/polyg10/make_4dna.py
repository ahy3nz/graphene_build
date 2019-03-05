import shutil
import os
import numpy as np
from numpy.random import randint

import script_utils
import operations
input_dir = '/raid6/homes/ahy3nz/Trajectories/charmm/polyg10/inputs/4dna'
base_dir = '/raid6/homes/ahy3nz/Trajectories/charmm/polyg10'
n_nodes = 2
jobid = [None] * n_nodes
def main():
    # For 0 dna systems
    trials = ['a' ,'b', 'c']
    constants = [125, 250]
    #constants = [50, 500]
    i = 0
    for constant in constants:
        for trial in trials:
            change_dir = os.path.join(
                    base_dir, '4dna', 'k{}_{}'.format(constant,trial))
            if not os.path.isdir(change_dir):
                os.makedirs(change_dir)
            os.chdir(change_dir)
            shutil.copyfile(os.path.join(input_dir, 'em.mdp'), './em.mdp')
            shutil.copyfile(os.path.join(input_dir, 'bilayer-graphene-dna-water-ions.gro'), './bilayer-graphene-dna-water-ions.gro')
            shutil.copyfile(os.path.join(input_dir, 'index.ndx'), './index.ndx')
            shutil.copyfile(os.path.join(input_dir, 'topol.top'), './topol.top')
            with open('angled_insertion.mdp', 'w') as f:
                _write_mdp(f, k=constant, seed=randint(100))
            with open('submit.pbs', 'w') as f:
                lines = 'cd {} \n'.format(change_dir)
                lines += 'gmx grompp -f em.mdp -c bilayer-graphene-dna-water-ions.gro -p topol.top -o em -maxwarn 1\n'
                lines += 'gmx mdrun -deffnm em \n'
                lines += 'gmx grompp -f angled_insertion.mdp -c em.gro -p topol.top -maxwarn 1 -n index.ndx -o pull \n'
                lines += 'gmx mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm pull \n'
                script_utils.write_rahman_script(f, 
                        jobname='4dna_k{}_{}'.format(constant,trial),
                        body=lines)
            operations.submit_job('submit.pbs', jobid, n_nodes, i)
            i += 1

def _write_mdp(f,k=500, seed=1):
    f.write("""
integrator               = md
dt                       = 0.001
nsteps                   = 5000000 ; 5e6 steps * 1e-3 fs/step = 5e3 ps = 5 ns
comm-mode                = Linear
nstcomm                  = 1000
comm-grps                = System

nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 1000
nstcalcenergy            = 1000
nstenergy                = 1000
nstxout-compressed       = 1000

cutoff-scheme            = Verlet
nstlist                  = 40
ns_type                  = grid
pbc                      = xyz
periodic-molecules       = no
dispcorr                 = no

coulombtype              = pme
coulomb-modifier         = none
rcoulomb                 = 1.2
vdw-type                 = cut-off
vdw-modifier             = force-switch

; cut-off lengths       
rlist                    = 1.2
rvdw-switch              = 1.0
rvdw                     = 1.2
fourierspacing           = 0.12
pme-order                = 4
ewald-rtol               = 1e-04
ewald-geometry           = 3d
epsilon-surface          = 0

tcoupl                   = nose-hoover
nsttcouple               = -1
nh-chain-length          = 10
tc-grps                  = system
tau_t                    = 1.0 
ref_t                    = 300 

pcoupl                   = no
pcoupltype               = semiisotropic
nstpcouple               = -1
; tau-p                    = 10.0 10.0
compressibility          = 4.5e-5 4.5e-5 
ref-p                    = 1.01325 1.01325

gen_vel                  = yes
gen_temp                 = 300.0
gen_seed                 = {seed}
annealing = no

pull = yes
pull-nstxout = 100
pull-nstfout = 100
pull-ngroups = 1
pull-ncoords = 1

pull-group1-name = anchor2
pull-coord1-groups = 0 1
pull-coord1-type = umbrella
pull-coord1-geometry =  direction-periodic
pull-coord1-vec =  0.000 0.000 -1.000
pull-coord1-origin = 4.217 7.780 -0.22
pull-coord1-rate = 0
pull-coord1-k = {k}
pull-coord1-start = No
""".format(k= k, seed= seed))

if __name__ == "__main__":
    main()
