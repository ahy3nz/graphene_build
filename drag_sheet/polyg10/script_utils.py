
def write_rahman_script(f, jobname="JOBNAME", body=""):
    """ Write rahman pbs script

    Parameters
    ----------
    f : file to write to
    """
    f.write("""#!/bin/sh -l
#PBS -N {jobname}
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -q low
#PBS -m abe
#PBS -M ayang41@gmail.com
#module load gromacs/5.1.4
module load gromacs/2018.1
{body}
""".format(jobname=jobname, body=body))
    return f


def write_edison_script(f, jobname="JOBNAME", N="1", t="24:00:00", body=""):
    """ Write edison slurm script

    Parameters
    ----------
    f : file to write to
    """
    f.write("""#!/bin/bash -l
#SBATCH -p regular
#SBATCH -N {N}
#SBATCH -t {t}
#SBATCH -J {jobname}
#SBATCH -o {jobname}.out
#SBATCH -e {jobname}.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu
#SBATCH -L SCRATCH

module load cray-fftw/3.3.6.3
module load lammps/2017.03.31
export OMP_NUM_THREADS=1
{body}
""".format(N=N, t=t, jobname=jobname, body=body))
    return f

def write_accre_script(f, jobname="JOBNAME", N="1", ntasks="2", gpu="1", time="0-96:00:00", body=""):
    """ Write accre slurm script

    Parameters
    ----------
    f : file to write to
    """
    f.write("""#!/bin/bash
#SBATCH --nodes={N}
#SBATCH --account=mccabe_gpu
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --partition=pascal
#SBATCH --gres=gpu:{gpu}
#SBATCH --time={time}
#SBATCH --output={jobname}.o
#SBATCH --error={jobname}.e
#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name={jobname}
module restore hoomd-py3
export OMP_NUM_THREADS=1
{body}
""".format(N=N, ntasks=ntasks, gpu=gpu, time=time, jobname=jobname, body=body)) 
    return f
        


