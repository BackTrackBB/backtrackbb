# -*- coding: utf8 -*-
import os
import numpy as np
import itertools
from collections import defaultdict
from scipy.ndimage.interpolation import zoom
from tatka_modules.bp_types import Trigger,Pick
from tatka_modules.map_project import rect2latlon
from tatka_modules.grid_projection import sta_GRD_Proj
from tatka_modules.plot import bp_plot
from tatka_modules.mod_bp_TrigOrig_time import TrOrig_time
import cPickle as pickle


def run_BackProj(args):
    return _run_BackProj(*args)


def _run_BackProj(idd, config, st, st_CF, frequencies,
                  stations, coord_sta, GRD_sta,
                  coord_eq):

    t_bb = np.arange(config.start_t, config.end_t, config.t_overlap)
    t_b = t_bb[idd]
    t_e = t_b + config.time_lag

    time = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate
    time_env = np.arange(st_CF[0].stats.npts) / st_CF[0].stats.sampling_rate

    grid1 = GRD_sta.values()[0]
    nx, ny, nz = np.shape(grid1.array)
    stack_grid = np.zeros((nx,ny,nz),float)
    arrival_times = defaultdict(list)
    trig_time = defaultdict(list)
    bp_trig_time = defaultdict(list)

    nn = int(config.t_overlap)

    fs_env = st_CF[0].stats.sampling_rate
    sttime_env = st_CF[0].stats.starttime

    rec_start_time = st[0].stats.starttime

    k = 0
    Mtau = []
    for sta1, sta2 in itertools.combinations(stations, 2):
        x_sta1, y_sta1 = coord_sta[sta1]
        x_sta2, y_sta2 = coord_sta[sta2]

        distance = np.sqrt((x_sta1-x_sta2)**2 + (y_sta1-y_sta2)**2)

        if distance <= config.maxSTA_distance:
            k += 1
            if config.varWin_stationPair:
                tau_max = GRD_sta[sta1].get_value(GRD_sta[sta2].sta_x,
                                                  GRD_sta[sta2].sta_y,
                                                  GRD_sta[sta2].sta_z)
                Mtau.append(np.round(tau_max,1))
                t_e = t_b + np.round(tau_max,1)
            else:
                Mtau = 'none'
                tau_max = None

            stack_grid += sta_GRD_Proj(st_CF, GRD_sta, sta1, sta2, t_b, t_e, nn,
                                       fs_env, sttime_env, config,
                                       nx, ny, nz, arrival_times, tau_max)
    stack_grid /= k

    max_stack_grid = np.where(stack_grid == np.max(stack_grid))
    i_max = max_stack_grid[0][0]
    j_max = max_stack_grid[1][0]
    k_max = max_stack_grid[2][0]

    if config.cut_data:
        start_tw = config.cut_start + t_b
    else:
        start_tw = t_b
        config.cut_start = 0.

    if config.save_projGRID:
        print 'saving GRIDS with results'
        out_file = os.path.join('out_grid', 'out_' + str(t_b) + '.pkl')
        pickle.dump(stack_grid, open(out_file, "wb"))

    trigger = None
    if stack_grid[i_max, j_max, k_max] >= config.trigger:
        if config.max_subdivide is not None:
            #We "zoom" the grid in a region close to the maximum
            zoom_factor = float(config.max_subdivide)
            #TODO: slicing will not work if the maximum
            #is close to the grid edge!
            slice_factor = 4 #this is hardcoded, for now
            sf = int(slice_factor)
            slice_grid = stack_grid[i_max-sf:i_max+sf,
                                    j_max-sf:j_max+sf,
                                    k_max-sf:k_max+sf]
            zoom_slice_grid = zoom(slice_grid, zoom_factor)
            zoom_i_max, zoom_j_max, zoom_k_max =\
                    [a[0]/zoom_factor
                     for a in np.where(zoom_slice_grid == np.max(zoom_slice_grid))]
            i_max = zoom_i_max + i_max-sf
            j_max = zoom_j_max + j_max-sf
            k_max = zoom_k_max + k_max-sf
        xx_trig, yy_trig, zz_trig = grid1.get_xyz(i_max, j_max, k_max)

        trigger = Trigger()
        trigger.x, trigger.y, trigger.z = grid1.get_xyz(i_max, j_max, k_max)
        trigger.i, trigger.j, trigger.k = i_max, j_max, k_max
        trigger.max_grid = np.round(np.max(stack_grid), 4)
        trigger.beg_win = start_tw
        trigger.end_win = start_tw + config.time_lag
        trigger.center_win = start_tw + config.time_lag/2.
        trigger.lat, trigger.lon =\
                rect2latlon(trigger.x, trigger.y)

    ##------------------Origin time calculation------------------------------------------------------
        bp_origin_time, bp_trig_time =\
                TrOrig_time(config, stations, GRD_sta, xx_trig, yy_trig, zz_trig,
                            rec_start_time, arrival_times, trig_time)
        pick = Pick()
        trigger.origin_time = bp_origin_time
        pick.station = stations
        pick.arrival_type = config.wave_type
        if bp_origin_time:
            trigger.eventid = bp_origin_time.strftime("%Y%m%d_%H%M")+'A'
            pick.eventid = trigger.eventid
            for sta in stations:
                pick.theor_time.append(bp_trig_time[sta][1])
                pick.pick_time.append(bp_trig_time[sta][4])
                
        trigger.add_picks(pick)
    ##-----------------------------------------------------------------------------------------------
        print trigger

    ## Plotting------------------------------------------------------------------
    n1 = 0
    n22 = len(frequencies) - 1
    fq_str = str(np.round(frequencies[n1])) + '_' + str(np.round(frequencies[n22]))
    datestr = st[0].stats.starttime.strftime('%y%m%d%H')

    if config.plot_results == 'True' or\
            (config.plot_results == 'trigger_only' and trigger is not None):
        bp_plot(config, grid1, stack_grid,
                coord_eq, t_b, t_e, datestr, fq_str,
                coord_sta, st, stations, st_CF,
                time, time_env,
                frequencies, n1, n22, trigger, arrival_times, bp_trig_time, Mtau)

    return trigger
