import numpy as np
import random
import os
import subprocess
import mdtraj
import json
import script_utils
from scipy.spatial import distance
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import plot_ay
plot_ay.setDefaults()

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

######################
### RWMD functions ###
######################

def write_rwmd_lines(gro='npt_500ps.gro', top='compound.top', 
                    tc_groups=['non-water','water'], cooling_rate=1000, 
                    t_pairs=[(305, 385), (305,385)]):
    """ Write submission lines

    Parameters
    ---------
    gro : str
    top : str
    tc_groups : list of str
    cooling_rate : float, default 1000
        RWMD ceiling drops by 1 K every `cooling_rate` ps
    t_pairs: n_groups x 2
       (t_min_group, t_max_group)
        
    Returns
    -------
    lines : sttr
    
    """
    
    n_cooling = rwmd.write_rwmd_files(tc_groups=tc_groups, gro=gro, top=top, 
                                     cooling_rate=cooling_rate, t_pairs=t_pairs)
    lines = """
gmx2018 grompp -f heating_phase.mdp -c {gro} -p {top} -o heating_phase &> heating_phase.out
gmx2018 mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm heating_phase -cpi heating_phase_prev.cpt -append

gmx2018 grompp -f cooling_phase0.mdp -c heating_phase.gro -p {top} -t heating_phase_prev.cpt -o cooling_phase0
gmx2018 mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm cooling_phase0 -cpi cooling_phase0_prev.cpt -append &> cooling_phase0.out


for ((i=1; i<={n_cooling} ; i++))
do
    gmx2018 grompp -f cooling_phase${{i}}.mdp -c cooling_phase$((${{i}}-1)).gro -p {top} -t cooling_phase$((${{i}}-1))_prev.cpt -o cooling_phase${{i}}
    gmx2018 mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm cooling_phase${{i}} -cpi cooling_phase${{i}}_prev.cpt -append &> cooling_phase${{i}}.out
done

gmx2018 grompp -f npt.mdp -c cooling_phase$((${{i}}-1)).gro -p {top} -o npt > npt_grompp.out
""".format(**locals())
    return lines

def analysis_routine(trajfile, grofile, pdbfile):

    import json
    from collections import OrderedDict
    import bilayer_analysis_functions
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    traj = mdtraj.load(trajfile, top=grofile)
    traj_pdb = mdtraj.load(trajfile, top=pdbfile)
    topol = traj.topology

    # Compute system information
    lipid_tails, headgroup_dict = bilayer_analysis_functions.identify_groups(traj, 
            forcefield='charmm36')
    n_lipid = len([res for res in traj.topology.residues if not res.is_water])
    n_lipid_tails = len(lipid_tails.keys())
    n_tails_per_lipid = n_lipid_tails/n_lipid



    # Vectorized Calculations start here
    apl_avg, apl_std, apl_list = bilayer_analysis_functions.calc_APL(traj,n_lipid, blocked=True)
    np.savetxt('apl.dat', apl_list)

    angle_avg, angle_std, angle_list = bilayer_analysis_functions.calc_tilt_angle(traj, topol, lipid_tails, blocked=True)
    np.savetxt('angle.dat', angle_list)

    apt_avg, apt_std, apt_list = bilayer_analysis_functions.calc_APT(traj, apl_list, angle_list, n_tails_per_lipid, 
            blocked=True)
    np.savetxt('apt.dat', apt_list)

    s2_ave, s2_std, s2_list = bilayer_analysis_functions.calc_nematic_order(traj, blocked=True)
    np.savetxt('s2.dat', s2_list)

    headgroup_distance_dict = bilayer_analysis_functions.compute_headgroup_distances(traj, topol, headgroup_dict, blocked=True)
    Hpp_ave, Hpp_std, Hpp_list = bilayer_analysis_functions.calc_bilayer_height(traj, headgroup_distance_dict, blocked=True, anchor='DSPC')
    np.savetxt('height.dat', Hpp_list)

    offset_dict = bilayer_analysis_functions.calc_offsets(traj, headgroup_distance_dict, blocked=True, anchor='DSPC')

    d_a, d_t, d_b, bins, interdig_list,interdig_avg, interdig_std = \
        bilayer_analysis_functions.calc_density_profile(traj, topol, 
                                                        blocked=True)
    np.savetxt('idig.dat', interdig_list)
    ##print('Calculating hydrogen bonds...')
    ##hbond_matrix_avg, hbond_matrix_std, hbond_matrix_list, labelmap = bilayer_analysis_functions.calc_hbonds(traj, traj_pdb, topol, lipid_dict, headgroup_dict)
    #
    # Printing properties
    outpdf = PdfPages(('bilayeranalysis.pdf'))
    datafile = OrderedDict()
    datafile['trajectory'] = trajfile
    datafile['structure'] = grofile
    datafile['n_frames'] = traj.n_frames
    datafile['lipids'] = n_lipid
    datafile['tails'] = n_lipid_tails
    datafile['APL'] = OrderedDict()
    datafile['APL']['unit'] = str(apl_avg.unit)
    datafile['APL']['mean'] = float(apl_avg._value)
    datafile['APL']['std'] = float(apl_std._value)
    datafile['APT'] = OrderedDict()
    datafile['APT']['unit'] = str(apt_avg.unit)
    datafile['APT']['mean'] = float(apt_avg._value)
    datafile['APT']['std'] = float(apt_std._value)
    datafile['Bilayer Height'] = OrderedDict()
    datafile['Bilayer Height']['unit'] = str(Hpp_ave.unit)
    datafile['Bilayer Height']['mean'] = float(Hpp_ave._value)
    datafile['Bilayer Height']['std'] = float(Hpp_std._value)
    datafile['Tilt Angle'] = OrderedDict()
    datafile['Tilt Angle']['unit'] = str(angle_avg.unit)
    datafile['Tilt Angle']['Bilayer'] = OrderedDict()
    datafile['Tilt Angle']['Bilayer']['mean'] = float(angle_avg._value)
    datafile['Tilt Angle']['Bilayer']['std'] = float(angle_std._value)
    datafile['S2'] = OrderedDict()
    datafile['S2']['mean'] = s2_ave
    datafile['S2']['std'] = s2_std
    datafile['Interdigitation'] = OrderedDict()
    datafile['Interdigitation']['unit'] = str(interdig_avg.unit)
    datafile['Interdigitation']['mean'] = float(interdig_avg._value)
    datafile['Interdigitation']['std'] = float(interdig_std._value)

    datafile['Offset'] = OrderedDict()
    for key in offset_dict.keys():
        datafile['Offset']['unit'] = str(offset_dict[key][0].unit)
        datafile['Offset'][key] = OrderedDict()
        datafile['Offset'][key]['mean'] = float(offset_dict[key][0]._value )
        datafile['Offset'][key]['std'] = float(offset_dict[key][1]._value )
        #datafile['Offset (A)'][key] = [str(offset_dict[key][0]), str(offset_dict[key][1])]

    datafile['Tilt Angle']['Leaflet 1'] = OrderedDict()
    datafile['Tilt Angle']['Leaflet 1']['mean'] = float(np.mean(angle_list[:, 
                                            0 :int(np.floor(n_lipid_tails/2))])._value)
    datafile['Tilt Angle']['Leaflet 1']['std'] = float(np.std(angle_list[:, 
                                                0 :int(np.floor(n_lipid_tails/2))])._value)

    datafile['Tilt Angle']['Leaflet 2'] = OrderedDict()
    datafile['Tilt Angle']['Leaflet 2']['mean'] = float(np.mean(angle_list[:, 
                                                    int(np.floor(n_lipid_tails/2)):])._value)
    datafile['Tilt Angle']['Leaflet 2']['std'] = float(np.std(angle_list[:, 
                                                int(np.floor(n_lipid_tails/2)):])._value)
    #for row_label in labelmap.keys():
    #    for col_label in labelmap.keys():
    #        row_index = labelmap[row_label]
    #        col_index = labelmap[col_label]
    #        hbond_avg = hbond_matrix_avg[row_index, col_index]
    #        hbond_std = hbond_matrix_std[row_index, col_index]
    #        outfile.write('{:<20s}: {} ({})\n'.format(str(row_label+"-"+ col_label), hbond_avg, hbond_std))


    # Plotting

    fig1 = plt.figure(1)
    plt.subplot(3,2,1)
    plt.plot(apl_list)
    plt.title('APL')

    plt.subplot(3,2,2)
    plt.plot(np.mean(angle_list, axis=1))
    plt.title('Tilt Angle ($^o$)')

    plt.subplot(3,2,3)
    plt.plot(np.mean(apt_list,axis=1))
    plt.title('APT')

    plt.subplot(3,2,4)
    plt.plot(Hpp_list)
    plt.title('H$_{PP}$')

    plt.subplot(3,2,5)
    plt.plot(s2_list)
    plt.title('S2')

    plt.subplot(3,2,6)
    plt.plot(interdig_list)
    plt.title('Interdigitation (A)')

    plt.tight_layout()
    outpdf.savefig(fig1)
    plt.close()

    density_profile_top_avg = np.mean(d_t, axis = 0)
    density_profile_bot_avg = np.mean(d_b, axis = 0)
    density_profile_avg  = np.mean(d_a, axis=0)
    #
    #
    fig2 = plt.figure(2)
    plt.subplot(2,1,1)
    plt.plot(bins,density_profile_avg)
    plt.xlabel('Depth (nm)')
    plt.title('Density Profile (kg m$^{-3}$)')


    plt.subplot(2,1,2)

    #plt.plot(bins,density_profile_bot_avg)
    #plt.plot(bins,density_profile_top_avg)

    plt.hist(np.mean(angle_list[:, 0 : int(np.floor(n_lipid_tails/2))], axis=0)._value, 
                    bins=50,  
                    alpha=0.5, facecolor='blue', normed=True)
    plt.hist(np.mean(angle_list[:, int(np.floor(n_lipid_tails/2)) : ], 
                    axis=0)._value, bins=50,  
                    alpha=0.5, facecolor='red', normed = True)
    plt.title('Angle Distribution by Leaflet')
    plt.xlabel('Angle ($^o$)')

    plt.tight_layout()
    outpdf.savefig(fig2)
    plt.close()
    outpdf.close()
    with open('data.txt', 'w') as f:
        json.dump(datafile, f, indent=2)


def calc_penetration(traj, bilayer_center=3.26, frame=-1):
    relevant_atoms = [a.index for res in traj.topology.residues 
                                if 'GRA_5' in res.name 
                                for a in res.atoms]
    depth = np.mean(traj.xyz[frame,relevant_atoms,2])
    fluc = np.mean(traj.xyz[-100:, relevant_atoms,2], axis=1) - bilayer_center
    fluc = np.std(fluc)
    d_from_center = depth - bilayer_center
    return d_from_center, fluc


    
def compute_work(chdir):
    """ For a given simulation compute work profiles"""
    os.chdir(chdir)
    print("Computing work in {}...".format(chdir))
    if os.path.isfile('pull.gro') and not os.path.isfile('work_profile.dat'):
        if not os.path.isfile('pull_nopbc.xtc'):
            p = subprocess.Popen('echo 0 | gmx trjconv -f pull.xtc -s pull.tpr -pbc mol -o pull_nopbc.xtc', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()
        traj = mdtraj.load('pull_nopbc.xtc', top='em.gro')
        box_dims = traj.unitcell_lengths[0]
        pullf = np.loadtxt('pull_pullf.xvg', comments=['@', '#'])
        raw_force = pullf
        
        
        # Reduce amount of data in pullf to match len of traj
        block_size = 10
        remainder = pullf.shape[0] % block_size
        pullf = pullf[0:-remainder,1]
        n_blocks = int(pullf.shape[0]/block_size)
        pullf = pullf.reshape(n_blocks, block_size)
        pullf = np.mean(pullf, axis=1)
        
        graphene_atom = int(open('index.ndx', 'r').readlines()[-1]) - 1
        integrand = []
        dist_from_center = []
        bilayer_center = 3.26
        past_5nm = False
        for i, force in enumerate(pullf):
            dx = distance.euclidean(traj.xyz[i, graphene_atom, :], 
                                    traj.xyz[i+1, graphene_atom, :])
            d_from_center, fluc = calc_penetration(traj, frame=i)
            if d_from_center <= 5.00:
                past_5nm = True
            dist_from_center.append( d_from_center )
            # Note the square box
            if dx > box_dims[0]:
                dx -= box_dims[0]

            if past_5nm:
                integrand.append(force*dx)
            else:
                integrand.append(0.0)
        
        work_profile = np.cumsum(integrand)
        np.savetxt("work_distance_profile.dat", 
                np.column_stack((dist_from_center, work_profile)))
        all_times = traj.time[1:]
        
        fig, ax = plt.subplots(1,1)
        ax.plot(raw_force[:,0], raw_force[:,1], color='r')
        ax.set_ylabel(r"Force ($\frac{kJ}{mol*nm}$)", color='r')
        for tick in ax.get_yticklabels():
            tick.set_color('r')
        
        ax2 = ax.twinx()
        ax2.plot(all_times, work_profile, color='b')
        np.savetxt('work_profile.dat', np.column_stack((all_times, work_profile)))
        ax2.set_ylabel(r"Work ($\frac{kJ}{mol}$)", color='b')
        for tick in ax2.get_yticklabels():
            tick.set_color('b')
        
        ax.set_xlabel("Time (ps)")
        fig.tight_layout()
        fig.savefig("profile.png", transparent=True)
        plt.close(fig)

def compute_penetration(chdir):
    """ Compute penetration single simulation """
    os.chdir(chdir)
    print("Computing penetration in {}...".format(chdir))
    if os.path.isfile("pull_nopbc.xtc") and not os.path.isfile('penetration.json'):
        traj = mdtraj.load('pull_nopbc.xtc', top ='em.gro')
        d_from_center, fluc = calc_penetration(traj)
        #all_dists.append(d_from_center)
        #all_flucs.append(fluc)
        #mean_dist = np.mean(all_dists)
        #std_dist = np.sqrt(np.sum([a**2 for a in all_fluc]))
        
        to_df = {
                'distance': float(d_from_center),
                'err': float(fluc)
                }

        json.dump(to_df, open('penetration.json', 'w'), indent=4)
        #df = pd.DataFrame(data=to_df, columns=['distance', 'err'])
        #df.to_csv('penetration.csv', columns=[ 'distance', 'err'])

