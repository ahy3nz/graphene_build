import pdb
import numpy as np
import mdtraj
from scipy.spatial.distance import euclidean

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import plot_ay
plot_ay.setDefaults()

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
l, = ax.plot(np.arange(traj.n_frames), np.mean(dopc_distances, axis=1),
        label='avg dopc distance')
ax.fill_between(np.arange(traj.n_frames), 
                np.mean(dopc_distances,axis=1) - np.std(dopc_distances,axis=1),
                np.mean(dopc_distances,axis=1) + np.std(dopc_distances,axis=1),
                alpha=0.4, color=l.get_color())

l, = ax.plot(np.arange(traj.n_frames), np.mean(chol_distances, axis=1),
        label='avg chol distance')
ax.fill_between(np.arange(traj.n_frames), 
                np.mean(chol_distances,axis=1) - np.std(chol_distances,axis=1),
                np.mean(chol_distances,axis=1) + np.std(chol_distances,axis=1),
                alpha=0.4, color=l.get_color())


ax.set_ylabel("Distance (nm)")
ax.set_xlabel("Frame")
plot_ay.tidyUp(fig, ax, gridArgs={}, legendArgs={}, tightLayoutArgs={})
fig.savefig('distances.png')
