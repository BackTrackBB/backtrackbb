# -*- coding: utf8 -*-
import os
import numpy as np
import itertools
from collections import defaultdict
from scipy.ndimage.interpolation import zoom
from bp_types import Trigger
from map_project import rect2latlon
from grid_projection import sta_GRD_Proj
from plot import bp_plot
from mod_bp_TrigOrig_time import TrOrig_time
from NLLGrid import NLLGrid


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

    # Create stack grid using a time grid as model
    gr = GRD_sta.values()[0]
    stack_grid = NLLGrid(nx=gr.nx, ny=gr.ny, nz=gr.nz,
                         x_orig=gr.x_orig, y_orig=gr.y_orig, z_orig=gr.z_orig,
                         dx=gr.dx, dy=gr.dy, dz=gr.dz)
    stack_grid.type = 'STACK'
    stack_grid.proj_name = gr.proj_name
    stack_grid.ellipsoid = gr.ellipsoid
    stack_grid.orig_lat = gr.orig_lat
    stack_grid.orig_lon = gr.orig_lon
    stack_grid.first_std_paral = gr.first_std_paral
    stack_grid.second_std_paral = gr.second_std_paral
    stack_grid.map_rot = gr.map_rot
    stack_grid.init_array()

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
                Mtau = None
                tau_max = None

            stack_grid.array += sta_GRD_Proj(st_CF, GRD_sta, sta1, sta2, t_b, t_e, nn,
                                       fs_env, sttime_env, config,
                                       arrival_times, tau_max)
    stack_grid.array /= k

    i_max, j_max, k_max = np.unravel_index(stack_grid.array.argmax(), stack_grid.array.shape)

    if config.cut_data:
        start_tw = config.cut_start + t_b
    else:
        start_tw = t_b
        config.cut_start = 0.

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
        xx_trig, yy_trig, zz_trig = stack_grid.get_xyz(i_max, j_max, k_max)

        trigger = Trigger()
        trigger.x, trigger.y, trigger.z = stack_grid.get_xyz(i_max, j_max, k_max)
        trigger.i, trigger.j, trigger.k = i_max, j_max, k_max
        trigger.max_grid = np.round(stack_grid.max(), 4)
        trigger.beg_win = start_tw
        trigger.end_win = start_tw + config.time_lag
        trigger.center_win = start_tw + config.time_lag/2.
        trigger.lat, trigger.lon =\
                rect2latlon(trigger.x, trigger.y)

    ##------------------Origin time calculation------------------------------------------------------
        bp_origin_time, bp_trig_time =\
                TrOrig_time(config, stations, GRD_sta, xx_trig, yy_trig, zz_trig,
                            rec_start_time, arrival_times, trig_time)

        trigger.origin_time = bp_origin_time
    ##-----------------------------------------------------------------------------------------------
        print trigger

    if config.save_projGRID == True or\
            (config.save_projGRID == 'trigger_only' and trigger is not None):
        print 'Saving projection grid to file.'
        basename = os.path.join(config.out_dir, 'out_t%05.1f' % t_b)
        stack_grid.write_hdr_file(basename)
        stack_grid.write_buf_file(basename)

    ## Plotting------------------------------------------------------------------
    n1 = 0
    n22 = len(frequencies) - 1
    fq_str = str(np.round(frequencies[n1])) + '_' + str(np.round(frequencies[n22]))
    datestr = st[0].stats.starttime.strftime('%y%m%d%H')
    bp_plot(config, stack_grid,
            coord_eq, t_b, t_e, datestr, fq_str,
            coord_sta, st, stations, st_CF,
            time, time_env,
            frequencies, n1, n22, trigger, arrival_times, bp_trig_time, Mtau)

    return trigger
