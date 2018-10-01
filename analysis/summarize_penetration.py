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
    if traj.n_frames > 1:
        frame = 1000
    else:
        frame = 0
    frame = -1
    depth = np.mean(traj.xyz[frame,relevant_atoms,2])
    d_from_center = depth - bilayer_center
    return d_from_center

def main():
    #dna_loading = ['alt6_dna']
    dna_loading = ['0_dna', '2_dna', '4_dna', 'alt6_dna']
    #constants = ['50']
    constants = ['50', '125']
    angles = ['0']
    #angles = ['0', '15', '30', '45']
    trials = ['a', 'b', 'c']
    curr_dir = os.getcwd()
    all_dicts = []
    for combo in itertools.product(dna_loading, constants, angles):
        all_dists = []
        for trial in trials:
            os.chdir(os.path.join(curr_dir, '../{0}/k{1}_{2}_{3}'.format(
                                                combo[0], combo[1], 
                                                combo[2], trial)))
            if os.path.isfile('pull.xtc') and os.path.isfile('em.gro'):
                #print(os.getcwd())
                try:
                    traj = mdtraj.load('pull.xtc', top ='em.gro')
                    print((os.getcwd(), traj.n_frames))
                    all_dists.append(calc_penetration(traj))
                except:
                    print((os.getcwd(), 'ERROR'))

        to_df = {'angle': combo[2],
                'constant': combo[1],
                'dna': combo[0],
                'distance': np.mean(all_dists),
                'err': np.std(all_dists)}
        all_dicts.append(to_df)

    os.chdir(curr_dir)
    df = pd.DataFrame(data=all_dicts)
    df.to_csv('penetration_5ns.csv', columns=['dna', 'constant', 'angle', 
                                        'distance', 'err'])

if __name__ == "__main__":
    main()
