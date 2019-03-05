import numpy as np
import pdb
import pandas as pd
import mdtraj
import itertools
import os

def calc_penetration(traj, bilayer_center=3.26):
    relevant_atoms = [a.index for res in traj.topology.residues 
                                if 'GRA_5' in res.name 
                                for a in res.atoms]
    frame = -1 
    depth = np.mean(traj.xyz[frame,relevant_atoms,2])
    d_from_center = depth - bilayer_center
    return d_from_center

def main():
    constants = ['50','500']
    trials = ['a', 'b', 'c']
    curr_dir = os.getcwd()
    all_dicts = []
    for constant in constants:
        all_dists = []
        for trial in trials:
            os.chdir(os.path.join(curr_dir, 'k{}_{}'.format(constant, trial))
            if os.path.isfile('pull_nopbc.xtc') and os.path.isfile('em.gro'):
                print(os.getcwd())
                traj = mdtraj.load('pull_nopbc.xtc', top ='em.gro')
                all_dists.append(calc_penetration(traj))

        to_df = {
                'constant': constant,
                'distance': np.mean(all_dists),
                'err': np.std(all_dists)
                'dna': '0dna'
                }
        all_dicts.append(to_df)

    os.chdir(curr_dir)
    df = pd.DataFrame(data=all_dicts)
    df.to_csv('penetration.csv', columns=['constant', 
                                        'distance', 'err', 'dna'])

if __name__ == "__main__":
    main()
