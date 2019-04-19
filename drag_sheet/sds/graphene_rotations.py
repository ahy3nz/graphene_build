import numpy as np
import sys
import mbuild as mb
import mdtraj
import pdb
import argparse

###########################
## This code is intended to transform the graphene + dna structure
## according to a specified angle of attack
## Additionally, the corresponding pulling MDP file is written
## based on this rotated-graphene sheet, such that the compound is pulled
## along that same vector/direction
############################

def align_sheet(cmpd_vector, ref_vector, cmpd):
    """ 
    Rotate sheet so east-west diagonal is parallel to X axis
    and sheet lies in XZ plane
    And sheet is now pointing down at the bilayer

    Parameters
    ---------
    cmpd_vector : 3-tuple
        Vector that represents graphene sheet diagonal
    ref_vector : 3-tuple
        Vector with which we are aligning
    cmpd : mb.Compound
    """
    theta = mb.coordinate_transform.angle(cmpd_vector, ref_vector)
    normal = np.cross(cmpd_vector, ref_vector)
    cmpd.rotate(theta, normal)

def calculate_pull_pts(point1, point2,  bilayer_center=3.26):
    """ Compute the pull vector and pulling reference points

    Parameters
    ---------
    point1 : 3-tuple
    point2 : 3-tuple
    init_distance : float (nm), default 2
        Distance (along the pull vector) between each graphene point
        and anchor point

    Returns
    ------
    pull_vec : 3-tuple
    anchor1 : 3-tuple
    anchor2 : 3-tuple

    Notes
    -----
    Anchors are stationary points within the bilayer that straddle the bilayer
    center.
    """
    ns_vector = point2 - point1
    ns_vector /= np.linalg.norm(ns_vector)

    # Identify the anchor points that properly
    # straddle the bilayer center while 
    # still being in line with the north-south vector of the graphene sheet
    z_gap = abs(point2[2] - point1[2])
    anchor1_z = bilayer_center + (z_gap/2)
    scale1_pull_vec = ( anchor1_z - point1[2] ) / ns_vector[2]
    anchor2_z = bilayer_center - (z_gap/2)
    scale2_pull_vec = ( anchor2_z - point2[2] ) / ns_vector[2]

    pull_vec1 = scale1_pull_vec * ns_vector
    pull_vec2 = scale2_pull_vec * ns_vector

    anchor1 = point1 + pull_vec1
    anchor2 = point2 + pull_vec2

    return ns_vector, anchor1, anchor2

def write_pull_mdp(pull_vec, anchor1, anchor2, filename='angled_insertion.mdp',
                    k=500):
    """ Write pulling mdp for this angle 

    Parameters
    ---------
    pull_vec : 3-tuple
    anchor1 : 3-tuple
    anchor2 : 3-tuple
    k : float, 
        Spring constant
    """

    with open('{}'.format(filename), 'w') as f:
        f.write("""integrator               = md
dt                       = 0.001
nsteps                   = 5000000 ; 5e6 steps * 1e-3 fs/step = 5e3 ps = 5 ns
;nsteps                   = 1000000 ; 1e6 steps * 1e-3 fs/step = 1e3 ps = 1 ns
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
gen_seed                 = 69
annealing = no

pull = yes
pull-nstxout = 100
pull-nstfout = 100
pull-ngroups = 1
pull-ncoords = 1
;pull-group1-name = anchor1
;pull-coord1-groups = 0 1
;pull-coord1-type = umbrella
;pull-coord1-geometry =  direction-periodic
;pull-coord1-vec = {:4.3f} {:4.3f} {:4.3f} 
;pull-coord1-origin = {:4.3f} {:4.3f} {:4.3f} 
;pull-coord1-rate = 0
;pull-coord1-k = 1000
;pull-coord1-start = No

pull-group1-name = anchor2
pull-coord1-groups = 0 1
pull-coord1-type = umbrella
pull-coord1-geometry =  direction-periodic
pull-coord1-vec =  {:4.3f} {:4.3f} {:4.3f}
pull-coord1-origin = {:4.3f} {:4.3f} {:4.3f}
pull-coord1-rate = 0
pull-coord1-k = {}
pull-coord1-start = No""".format(*pull_vec, *anchor1, *pull_vec, *anchor2, k))
    

def rotate_gnf_write_mdp(gro_file='composite.gro', angle_of_attack=0,
        force_constant=50, outgro='composite_rotated.gro'):
    """ Rotate GNF and write pulling MDP file"""
    cmpd = mb.load(gro_file)
    traj = mdtraj.load(gro_file)

    old_pos = cmpd.pos
    
    # Define vectors
    
    east = cmpd.xyz[80]
    west = cmpd.xyz[926]
    cmpd_vector = east - west
    ref_vector = [1,0,0]
    
    # Align the graphene sheet and reset its orientation
    align_sheet(cmpd_vector, ref_vector, cmpd)
    
    # Now spin to change angle of attack
    cmpd.rotate(np.deg2rad(angle_of_attack), [1,0,0])
    
    # Update coordinates and save grofile via mdtraj
    cmpd.translate_to(old_pos)
    print("Copying atoms over")
    #for i in range(traj.n_atoms):
        #pdb.set_trace()
    traj.xyz[0,:] = cmpd.xyz[:]
    traj.save(outgro)
    
    # Get pulling points and references for MDP file
    north = cmpd.xyz[0]
    south = cmpd.xyz[1006]
    pull_vec, anchor1, anchor2 = calculate_pull_pts(north, south, 
                                                    bilayer_center=3.26)
    
    # Write the mdp file 
    write_pull_mdp(pull_vec, anchor1, anchor2,k=force_constant)


