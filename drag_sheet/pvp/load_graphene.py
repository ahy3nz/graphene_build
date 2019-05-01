import os 
from itertools import groupby
import numpy as np

import mbuild as mb
from mbuild.packing import fill_box
import mdtraj

import graphene_rotations

PATH_TO_MOLECULES = '/raid6/homes/ahy3nz/Programs/graphene_build/single_molecules/'
def construct_system(dist_from_sheet=0.7, box_thickness=0.5, n_load_per_side=10,
        bilayer_center=[5.387, 5.347, 3.26]):
    """ Build graphene + PVP + bilayer + ions
    The idea is to pack a box below the sheet and below the sheet.
    The box we are filling is about 5nm x 5nm x 0.635nm. 

    Parameters
    ---------
    dist_from_sheet : float, default 0.7 nm
        When defining the box above/below the sheet, specify how far above/below
    box_thickness : float: default 0.5nm
        When defining the box to fill, the X and Y dimensions take from
        the graphene sheet, but the Z dimensions (thickness) is defined 
        in box_thickneess
    n_load_per_side : int, default 10
        Important for packing the boxes and specifying number of PVP 

    """

    header="""; Include forcefield parameters
#include "/raid6/homes/ahy3nz/Programs/graphene_build/force-field/ffmix.itp"
#include "/raid6/homes/ahy3nz/Programs/graphene_build/force-field/molecules.itp"
#include "/raid6/homes/ahy3nz/Programs/graphene_build/force-field/ions.itp
    """


    gra = mb.load(os.path.join(PATH_TO_MOLECULES, 'sheet.gro'))
    gra.name = 'GRA_5'

    pvp = mb.load(os.path.join(PATH_TO_MOLECULES, 'PVP.gro'))
    pvp.name = 'PVP'
    # The PVP  molecules are initially aligned with the Z-axis, so do a rotation
    pvp.spin(np.pi/2, [1,0,0]) 

    
    above_box = mb.Box(mins=np.min(gra.xyz, axis=0)+[0, 0, dist_from_sheet], 
                       maxs=np.max(gra.xyz, axis=0)+[0, 0, dist_from_sheet + box_thickness])
    fill_above = fill_box(pvp, 
            n_compounds=n_load_per_side, box=above_box,
            edge=0.1, overlap=0.1)

    below_box = mb.Box(mins=np.min(gra.xyz, axis=0)-[0, 0, dist_from_sheet + box_thickness], 
                       maxs=np.max(gra.xyz, axis=0)-[0, 0, dist_from_sheet])
    fill_below = fill_box(pvp, 
            n_compounds=n_load_per_side, box=below_box,
            edge=0.1, overlap=0.1)

    sys = mb.Compound()
    sys.add([gra, fill_above, fill_below])

    sys.spin(np.pi/2, [1,0,0]) # Do some rotations to get it corner first
    sys.spin(np.pi/4, [0,1,0])
    sys.translate_to(bilayer_center) # This is the bilayer center
    sys.translate([0,0,6]) # Shift the loaded GNF up to be above the bialyer

    # For the sake of writing, we will be using a lot of mdtraj functionality
    traj = sys.to_trajectory(residues=['GRA_5', 'PVP'])
    bilayer = mdtraj.load(os.path.join(PATH_TO_MOLECULES, 'bilayer.gro'))
    composite_traj = traj.stack(bilayer)
    composite_traj.unitcell_lengths = bilayer.unitcell_lengths
    composite_traj.unitcell_lengths[0,2] = 13
    composite_traj[0].save('composite.gro')
    write_gmx_topology('topol.top', composite_traj, header=header)

    # Use external module to perform rotations and write pulling MDP
    graphene_rotations.rotate_gnf_write_mdp(gro_file='composite.gro', 
            angle_of_attack=0,
            force_constant=50, outgro='composite_rotated.gro')



def write_gmx_topology(outfile, composite_traj, header=""):
    with open(outfile, 'w') as f:
        f.write("{}\n\n".format(header))
        f.write("[ system ]\n")
        f.write("Loaded-GNF-bilayer system\n\n")
        f.write("[ molecules ] \n")
        f.write("\n".join(["{:<8s} {}".format(name, sum(1 for _ in g))
                          for name, g in groupby([c.name for c in 
                              composite_traj.topology.residues])]))
        f.write("\n")

if __name__ == "__main__":
    construct_system()
