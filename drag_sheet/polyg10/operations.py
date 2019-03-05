import numpy as np
import random
import os
import subprocess
import mdtraj

def modify_top(top='compound.top',
               include="#include \"/raid6/homes/ahy3nz/Programs/McCabeGroup/atomistic/forcefield.itp\""):
    """ Modify topology file include statement

    Parameters
    ---------
    top : str
        filename of gmx topology
    include : str
        new include statement

    Notes'
    ------
    Top file is overwritten
    Assumes the include statement is in the first line

    Useful includes to remember:
    rahman: /raid6/homes/ahy3nz/Programs/McCabeGroup/atomistic/forcefield.itp
    edison: /global/homes/a/ahy3nz/Programs/McCabeGroup/atomistic/forcefield.itp
    accre: /home/yangah/Programs/McCabeGroup/atomistic/forcefield.itp

    """
    toplines = open(top,'r').readlines()
    toplines[0] = include + "\n"
    with open(top,'w') as f:
        for line in toplines:
            f.write(line)

def submit_job(script, jobid, n_nodes, i):
    """ Submit job to cluster 

    Parameters
    ---------
    script : str
        name of submission script
    n_nodes : int
        number of nodes we are using, equivalently the number
        of jobs running at any one time with others on hold
    jobid : list 
        Constains the last submited job for each node
    i : int
        Index/iteration number of the job we are submitting

    Returns
    -------
    jobid : list
        updated list corresponding to jobids last submitted


    """
    if 'pbs' in script:
        if not jobid[i % n_nodes]:
            p = subprocess.Popen('qsub {}'.format(script), shell=True, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)
            p.wait()
        elif jobid[i % n_nodes]:
            p = subprocess.Popen('qsub -W depend=afterany:{0} {1}'.format(\
                                 jobid[i % n_nodes], script), shell=True, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)
            p.wait()
        outs, errs = p.communicate()
        jobid[i % n_nodes] = outs.decode().strip()
    elif "sbatch" in script:
        if not jobid[i % n_nodes]:
            p = subprocess.Popen('sbatch {}'.format(script), shell=True, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)
            p.wait()
        elif jobid[i % n_nodes]:
            p = subprocess.Popen('sbatch -d afterany:{0} {1}'.format(\
                                 jobid[i % n_nodes], script), shell=True, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)
            p.wait()
        outs, errs = p.communicate()
        jobid[i % n_nodes] = outs.decode().strip().split()[-1]
    return jobid

###############################
### Equilibration functions ###
###############################

def write_eq_lines(gro='compound.gro', top='compound.top'):
    """ Write EM, NVT, NPT lines for equilibration """
    lines = """
gmx grompp -f em.mdp -c {gro} -p {top} -o em -maxwarn 2 &> em_grompp.log
gmx mdrun -deffnm em 

gmx grompp -f nvt.mdp -c em.gro -p {top} -o nvt -maxwarn 2 &> nvt_grompp.log
gmx mdrun -deffnm nvt

gmx grompp -f npt_500ps.mdp -c nvt.gro -p {top} -o npt_500ps -t nvt.cpt -maxwarn 2 &> npt_500ps_grompp.log
gmx mdrun -deffnm npt_500ps -ntomp 8 -ntmpi 2 -gpu_id 01 """.format(**locals())
    return lines

###########################
## Production functions ###
###########################

def write_production_lines(filename='npt'):
    """ Write NPT production """
    lines = "gmx mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm {filename} -cpi {filename}_prev.cpt -append".format( **locals())
    return lines


