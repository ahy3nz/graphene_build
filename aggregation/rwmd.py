import numpy as np
import tex_sampling.pyrex as pyrex


######################
### RWMD functions ###
######################
### A random-walk through temperature space, attempting temperature swaps with
### adjacent temperature windows every so often. Includes functions for 
### a ramp-down cycle when the RWMD temperature ceiling decreases

def write_rwmd_files(tc_groups=['non-water','water'], gro='npt_500ps.gro', 
                    top='compound.top', 
                    cooling_rate=1000, t_pairs=[[305, 385], [305,385]]):
    """ Write files relevant for RWMD

    Parameters
    ---------
    components : []
        List of each component (str) in system 
    gro : str
    top : str
    cooling_rate : float, default 1000
        RWMD ceiling drops by 1 K every `cooling_rate` ps
    t_pairs: n_groups x 2
       (t_min_group, t_max_group)


    Returns
    ------
    n_cooling : int
        last # of cooling cycle
        n_cooling 4 means the last index is 4 (so there are actually 5 cooling phases)


    Notes
    -----
    Will require a recompiled verion of gromacs, modified for longer string buffers
    """
    # First 25 ns is just the maximum range of heating
    t_pairs = np.asarray(t_pairs)
    assert np.shape(t_pairs)[0] == len(tc_groups), "tc_groups needs to match t_pairs"

    current_T, t_pairs, sim_time = tex_heating(t_pairs=t_pairs, tc_groups=tc_groups,
                                                dT=10, current_T=305, sim_time=0,
                                                interval_duration=5, timestep=0.002,
                                                heating_duration=30000)

    n_cooling = tex_cooling(t_pairs=t_pairs, current_T=current_T, 
                sim_time=sim_time,
                cooling_duration=30000, interval_duration=5, cooling_rate=1000,
                timestep=0.002, tc_groups=tc_groups)

    #return _write_rahman_rwmd(gro=gro, top=top, n_cooling=n_cooling)
    return n_cooling


def tex_heating(t_pairs=[(305,385), (305,385)], dT=10, current_T=305, 
                sim_time=0, heating_duration=30000, interval_duration=5,
                timestep=0.002, tc_groups=['non-water', 'water']):
    """ Temperature exchange heating phase

    Parameters
    ---------
    t_pairs: n_groups x 2
       (t_min_group, t_max_group)
    dT : int, optional, default 10 (K)
        temperature spacing
    current_T: float, default 305 (K)
        starting temperature
    sim_time : float, default 0 (ps)
        starting simulation time
    heating_duration : float, default 30000 (ps)
        Duration of heating phase
    interval_duration: float, default 5 (ps)
        Duration of sampling at a particular temperature window before attempting
        a temperature exchange

    Returns
    -------
    annealing_times: list
        List of times for gmx annealing
    temps: # temperature coupling groups x # temperature points
        List of temperatures for gmx annealing, each row corresponds to each 
        temperature coupling group
    t_pairs: # tempreature coupling groups x 2
        (t_min_group, t_max_group)
    sim_time:  float (ps)
        Simulation time at the end of the heating phase


    Notes
    -----
    Observe that the hard-coded timestep is 0.002 ps per timestep
    Most loops and counters are based on timesteps and not explicit times
    """
    t_pairs = np.asarray(t_pairs)
    n_groups = np.shape(t_pairs)[0]
    timestep = 0.002 # 2fs / timestep

    # The full heating phase
    heating_steps = int(heating_duration / timestep) # steps
    interval_steps = int(interval_duration / timestep)  #steps
    for step in range(0,heating_steps):
        sim_time += timestep
        # Initialization
        if step == 0:
            freq = pyrex.init_freq_dict(t_pairs[0,0], t_pairs[0,1], dT, T_init=current_T)
            temps = np.zeros((n_groups, int(heating_steps/interval_steps)) )
            temps[:,0] = current_T
            annealing_times = np.zeros(int(heating_steps/interval_steps), dtype=int)
            annealing_times[0] = int(sim_time) 
            _adjust_t_groups(t_pairs, temps, current_T, step, interval_steps,
                        thermostat_style='plateau')
        # Attempt to change temperature
        if step % interval_steps == 0:
            current_T = pyrex.choose_next_T(freq, current_T, dT)
            temps[0, int(step/interval_steps)] = current_T
            annealing_times[int(step/interval_steps)] = int(sim_time)
            _adjust_t_groups(t_pairs, temps, current_T, step, interval_steps,
                        thermostat_style='plateau')

    with open("heating_phase.mdp", 'w')as f:
        _write_body(f, ref_temp=t_pairs[0,0], n_steps=heating_steps, 
                    tc_groups=tc_groups)
        _write_annealing_lines(f, annealing_times, temps)


    return current_T, t_pairs, sim_time

def tex_cooling(t_pairs=[(305,385), (305,385)], dT=10, current_T=305, 
                sim_time=30000, cooling_duration=30000, interval_duration=5,
                cooling_rate=1000, timestep=0.002,
                thermostat_style='plateau', tc_groups=['non-water','water']):
    """ Temperature exchange heating phase

    Parameters
    ---------
    t_pairs: n_groups x 2
       (t_min_group, t_max_group)
    dT : int, optional, default 10 (K)
        temperature spacing
    current_T: float, default 305 (K)
        starting temperature
    sim_time : float, default 0 (ps)
        starting simulation time
    cooling_duration : float, default 30000 (ps)
        Duration of a single cooling simulation
    interval_duration: float, default 5 (ps)
        Duration of sampling at a particular temperature window before attempting
        a temperature exchange
    timestep : float, default 0.002 (ps/step)
    cooling_rate : float, default 1000
        Every 1000 ps, reduce temperature ceiling by 1 K
    thermostat_style: str, default='plateau'
       'plateau' : If the primary group temperature exceeds the maximum temperature
            for other groups, set the other group's temperatures to their respective
            maxima 
        'proportional' : Temperature changes from the primary group are scaled
            for each temperature group


    Returns
    -------
    annealing_times: list
        List of times for gmx annealing
    temps: # temperature coupling groups x # temperature points
        List of temperatures for gmx annealing, each row corresponds to each 
        temperature coupling group
    t_pairs: # tempreature coupling groups x 2
        (t_min_group, t_max_group)
    sim_time:  float (ps)
        Simulation time at the end of the heating phase


    Notes
    -----
    Observe that the hard-coded timestep is 0.002 ps per timestep
    Most loops and counters are based on timesteps and not explicit times
    cooling_rate is extrapolated to eliminate one of the temperature windows
    """

    # Here's the cooling phase, but implemented as multiple cooling simulations
    # Extrapolate cooling rate and convert to frames
    t_pairs = np.asarray(t_pairs)
    n_groups = np.shape(t_pairs)[0]
    cooling_rate_steps = cooling_rate * dT / timestep 
    total_cooling_time = (t_pairs[0,1] - t_pairs[0,0]) * cooling_rate 

    n_cooling_phases = int(np.ceil(total_cooling_time/cooling_duration))
    cooling_steps = int(cooling_duration / timestep)#steps
    interval_steps = int(interval_duration / timestep)#steps

    for i in range(n_cooling_phases):
        t_init = sim_time
        for step in range(0, cooling_steps):
            sim_time += timestep
            # Initialization
            if step == 0:
                temps = np.zeros((n_groups, int(cooling_steps/interval_steps)))
                temps[:,0] = current_T
                annealing_times = np.zeros(int(cooling_steps/interval_steps))
                annealing_times[0] = int(sim_time)
                freq = pyrex.init_freq_dict(t_pairs[0,0], t_pairs[0,1], 
                                            dT, T_init=current_T)
            # Attempt to drop the RWMD ceiling temperature
            if step % cooling_rate_steps  == 0 and t_pairs[0,1] > t_pairs[0,0] + dT:
                for t_index in range(t_pairs.shape[0]):
                    if t_pairs[t_index, 1] > t_pairs[t_index,0] + dT:
                        t_pairs[t_index,1] -= dT
                if t_pairs[0, 1] < current_T:
                    current_T = t_pairs[0, 1]
                freq = pyrex.init_freq_dict(t_pairs[0,0], t_pairs[0,1], 
                                            dT, T_init=current_T)
            # Attempt to change temperature
            if step % interval_steps == 0:
                current_T = pyrex.choose_next_T(freq, current_T, dT)
                temps[0, int(step/interval_steps)] = current_T
                annealing_times[int(step/interval_steps)] = int(sim_time)
                _adjust_t_groups(t_pairs, temps, current_T, step, interval_steps,
                             thermostat_style=thermostat_style)
        with open("cooling_phase{}.mdp".format(i),'w') as f:
            _write_body(f, ref_temp=temps[0,0], n_steps=cooling_steps, 
                        t_init=t_init, tc_groups=tc_groups)
            _write_annealing_lines(f, annealing_times, temps)

    return i


def _write_body(f, ref_temp=305, n_steps=50000000,t_init=0,tc_groups=None):
    tc_grps_string = " ".join([thing for thing in tc_groups])
    tau_t_string = " ".join(["1.0" for _ in range(len(tc_groups))])
    ref_t_string = " ".join(["{:5.0f}".format(ref_temp) for _ in range(len(tc_groups))])
    annealing_string = " ".join(["single" for _ in range(len(tc_groups))])
    f.write(""" title                       = RWMD
; Run parameters
integrator                  = md
nsteps                      = {n_steps}     ; 2fs/step * 5e7 steps = 5e7 fs = 50 ns
dt                          = 0.002
tinit                       = {t_init:.0f}

; Output control
nstxout                     = 0             ; Don't save coordinates 
nstvout                     = 0             ; Don't save velocities
nstenergy                   = 10000
nstlog                      = 10000
nstxtcout                   = 10000

;bond parameters
continuation                = yes
constraint_algorithm        = lincs
constraints                 = h-bonds
;constraints                 = none
lincs_iter                  = 1
lincs_order                 = 4

; Neighbor searching
cutoff-scheme               = Verlet
ns_type                     = grid
nstlist                     = 10
rcoulomb                    = 1.4
rvdw                        = 1.4

;Electrostatics
coulombtype                 = PME
pme_order                   = 4
fourierspacing              = 0.16

; Temperature coupling
tcoupl                      = nose-hoover
tc-grps                     = {tc_grps_string}    
tau_t                       = {tau_t_string}
ref_t                       = {ref_t_string}

;Pressure coupling
pcoupl                      = Parrinello-Rahman
pcoupltype                  = semiisotropic
tau_p                       = 5.0           ; ps
ref_p                       = 1.0 1.0          ; bar   
compressibility             = 4.5e-5 4.5e-5
refcoord_scaling            = com

;PBC
pbc                         = xyz

;Dispersion correction
DispCorr                    = EnerPres

;Velocity generation
gen_vel                     = no

;Simulated annealing
annealing                   = {annealing_string}
""".format(**locals()))


def _write_annealing_lines(f, times, all_temps):
    if len(times) == np.shape(all_temps)[1]:
        temp_string = ""
        small_time_string = "".join(["{} ".format(thing) for thing in times])
        for i in range(all_temps.shape[0]):
            temp_string += "".join(["{:3.0f} ".format(thing) for thing in all_temps[i,:]])


        n_points = "".join(["{} ".format(len(times)) for _ in range(all_temps.shape[0])])
        time_string = " ".join([small_time_string for _ in range(all_temps.shape[0])])


        f.write("annealing-npoints\t={0}\n".format(n_points))
        f.write("annealing-time\t={0}\n".format(time_string))
        f.write("annealing-temp\t={0}\n".format(temp_string))
    else:
        sys.exit("Error, annealing lines do not have matching dimensions")


def _adjust_t_groups(t_pairs, temps, current_T, step, interval_steps,
                    thermostat_style='plateau'):
    """ Adjust the other temperature groups 

    t_pairs: n_groups x 2
       (t_min_group, t_max_group)
    temps : list
    current_T : float
    step : current timestep the loop is on
    interval_steps : in units of timestep, the number of steps until we changea
    thermostat_style: str, default='plateau'
       'plateau' : If the primary group temperature exceeds the maximum temperature
            for other groups, set the other group's temperatures to their respective
            maxima 
        'proportional' : Temperature changes from the primary group are scaled
            for each temperature group
        
    """
    for t_group in range(1, np.shape(t_pairs)[0]):
        if thermostat_style == 'proportional':
            primary_temp_change = current_T - t_pairs[0][0]
            scaling_factor = ((t_pairs[t_group][1] - t_pairs[t_group][0]) /
                             (t_pairs[0][1] - t_pairs[0][0]))
            temps[t_group, int(step/interval_steps)] = (scaling_factor 
                                                        * primary_temp_change)
        elif thermostat_style == 'plateau':
            if current_T > t_pairs[t_group][1]:
                temps[t_group, int(step/interval_steps)] = t_pairs[t_group][1]
            else:
                temps[t_group, int(step/interval_steps)] = current_T

