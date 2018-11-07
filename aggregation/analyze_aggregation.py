import os
import subprocess
import numpy as np
import mdtraj
from scipy.spatial.distance import euclidean

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import plot_ay
plot_ay.setDefaults()

def main():
    curr_dir = os.getcwd()
    folders = ['trial1', 'trial2', 'trial3']
    for folder in folders:
        os.chdir(os.path.join(curr_dir, folder))
        cmd = 'echo 0 0 | gmx trjconv -f npt.xtc -s npt.tpr -pbc mol -o npt_nopbc.xtc'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        p.wait()
        analyze_time_aggregation(xtc='npt_nopbc.xtc', gro='npt.gro')

def analyze_time_aggregation(xtc='npt_nopbc.xtc', gro='npt.gro'):
    traj = mdtraj.load('npt_nopbc.xtc', top='npt.gro')
    gra_atoms = [a.index for r in traj.topology.residues if 'GRA_5' in r.name for
                    a in r.atoms]
    dopc_residues = [r for r in traj.topology.residues if 'DOPC' in r.name]
    chol_residues = [r for r in traj.topology.residues if 'CHL1' in r.name]
    
    
    gra_centers = mdtraj.compute_center_of_mass(traj.atom_slice(gra_atoms))
    dopc_centers = np.zeros((traj.n_frames, len(dopc_residues), 3))
    chol_centers = np.zeros((traj.n_frames, len(chol_residues), 3))
    
    dopc_distances = np.zeros((traj.n_frames, len(dopc_residues)))
    chol_distances = np.zeros((traj.n_frames, len(chol_residues)))
    
    frame_step = 5 # one frame is 5 ps
    
    for i, dopc in enumerate(dopc_residues):
        dopc_atoms = [a.index for a in dopc.atoms]
        dopc_centers[:, i, :] = mdtraj.compute_center_of_mass(traj.atom_slice(dopc_atoms))
    for i, chol in enumerate(chol_residues):
        chol_atoms = [a.index for a in chol.atoms]
        chol_centers[:, i, :] = mdtraj.compute_center_of_mass(traj.atom_slice(chol_atoms))
    
    for i, gra_center in enumerate(gra_centers):
        for j in range(len(dopc_residues)):
            chol_distances[i,j] = euclidean(gra_center, chol_centers[i,j,:])
        for j in range(len(chol_residues)):
            dopc_distances[i,j] = euclidean(gra_center, dopc_centers[i,j,:])
    
    
    fig, ax = plt.subplots(1,1)
    l, = ax.plot(frame_step * np.arange(traj.n_frames), np.mean(dopc_distances, axis=1),
            label='Avg dopc distance')
    ax.fill_between(frame_step * np.arange(traj.n_frames), 
                    np.mean(dopc_distances,axis=1) - np.std(dopc_distances,axis=1),
                    np.mean(dopc_distances,axis=1) + np.std(dopc_distances,axis=1),
                    alpha=0.4, color=l.get_color())
    np.savetxt('dopc_distances.dat',
                                np.column_stack((frame_step*np.arange(traj.n_frames),
                                np.mean(dopc_distances,axis=1),
                                np.std(dopc_distances,axis=1))))
    
    l, = ax.plot(frame_step * np.arange(traj.n_frames), np.mean(chol_distances, axis=1),
            label='Avg chol distance')
    ax.fill_between(frame_step * np.arange(traj.n_frames), 
                    np.mean(chol_distances,axis=1) - np.std(chol_distances,axis=1),
                    np.mean(chol_distances,axis=1) + np.std(chol_distances,axis=1),
                    alpha=0.4, color=l.get_color())
    np.savetxt('chol_distances.dat',
                                np.column_stack((frame_step*np.arange(traj.n_frames),
                                np.mean(chol_distances,axis=1),
                                np.std(chol_distances,axis=1))))

    
    
    ax.set_ylabel("Distance (nm)")
    ax.set_xlabel("Time (ps)")
    plot_ay.tidyUp(fig, ax, gridArgs={}, legendArgs={}, tightLayoutArgs={})
    fig.savefig('Distances.png')

if __name__ == "__main__":
    main()
