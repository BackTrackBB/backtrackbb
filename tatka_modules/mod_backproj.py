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


def _run_BackProj(config, st, st_CF, t_begin, frequencies,
                  coord_sta, GRD_sta, coord_eq):

    t_end = t_begin + config.time_lag

    # Create stack grid using a time grid as model
    gr = GRD_sta.values()[0].values()[0]
    stack_grid = NLLGrid(nx=gr.nx, ny=gr.ny, nz=gr.nz,
                         x_orig=gr.x_orig, y_orig=gr.y_orig, z_orig=gr.z_orig,
                         dx=gr.dx, dy=gr.dy, dz=gr.dz)
    stack_grid.type = 'STACK'
    if gr.proj_name != 'NONE':
        stack_grid.proj_name = gr.proj_name
        stack_grid.ellipsoid = gr.ellipsoid
        stack_grid.orig_lat = gr.orig_lat
        stack_grid.orig_lon = gr.orig_lon
        stack_grid.first_std_paral = gr.first_std_paral
        stack_grid.second_std_paral = gr.second_std_paral
        stack_grid.map_rot = gr.map_rot
    stack_grid.init_array()

    arrival_times = defaultdict(lambda: defaultdict(list))

    sta_wave = [(sta, wave) for wave in config.wave_type for sta in config.stations]

    k = 0
    Mtau = []
    for sta_wave1, sta_wave2 in itertools.combinations(sta_wave, 2):
        sta1 = sta_wave1[0]
        sta2 = sta_wave2[0]
        wave1 = sta_wave1[1]
        wave2 = sta_wave2[1]

        x_sta1, y_sta1 = coord_sta[sta1]
        x_sta2, y_sta2 = coord_sta[sta2]
        distance = np.sqrt((x_sta1-x_sta2)**2 + (y_sta1-y_sta2)**2)
        if distance > config.maxSTA_distance:
            continue

        if config.varWin_stationPair:
            tau_max = GRD_sta[sta1][wave1].get_value(GRD_sta[sta2][wave1].sta_x,
                                                     GRD_sta[sta2][wave1].sta_y,
                                                     GRD_sta[sta2][wave1].sta_z)
            Mtau.append(np.round(tau_max,1))
            t_end = t_begin + np.round(tau_max,1)
        else:
            Mtau = None
            tau_max = None

        proj_grid, arrival1, arrival2 =\
                sta_GRD_Proj(config, st_CF, GRD_sta, sta1, sta2, wave1, wave2, t_begin, t_end, tau_max)
        stack_grid.array += proj_grid
        arrival_times[sta1][wave1].append(arrival1)
        arrival_times[sta2][wave2].append(arrival2)
        k += 1
    stack_grid.array /= k

    i_max, j_max, k_max = np.unravel_index(stack_grid.array.argmax(), stack_grid.array.shape)

    if config.cut_data:
        start_tw = config.cut_start + t_begin
    else:
        start_tw = t_begin
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
            # check if zoom_slice_grid is not empty.
            # TODO: why can it be empty?
            if zoom_slice_grid.size > 0:
                zoom_i_max, zoom_j_max, zoom_k_max =\
                     [a[0]/zoom_factor
                      for a in np.where(zoom_slice_grid == np.max(zoom_slice_grid))]
                i_max = zoom_i_max + i_max-sf
                j_max = zoom_j_max + j_max-sf
                k_max = zoom_k_max + k_max-sf

        trigger = Trigger()
        trigger.x, trigger.y, trigger.z = stack_grid.get_xyz(i_max, j_max, k_max)
        trigger.i, trigger.j, trigger.k = i_max, j_max, k_max
        trigger.max_grid = np.round(stack_grid.max(), 4)
        trigger.beg_win = start_tw
        trigger.end_win = start_tw + config.time_lag
        trigger.center_win = start_tw + config.time_lag/2.
        if stack_grid.proj_name != 'NONE':
            trigger.lat, trigger.lon =\
                    rect2latlon(trigger.x, trigger.y)

        # Compute origin time and theoretical arrivals
        TrOrig_time(config, GRD_sta, trigger, arrival_times)
        if trigger.origin_time is None:
            trigger = None

    if trigger is not None:
        print trigger

    if config.save_projGRID == True or\
            (config.save_projGRID == 'trigger_only' and trigger is not None):
        print 'Saving projection grid to file.'
        basename = os.path.join(config.out_dir, 'out_t%05.1f' % t_begin)
        stack_grid.write_hdr_file(basename)
        stack_grid.write_buf_file(basename)

    ## Plotting------------------------------------------------------------------
    n1 = 0
    n22 = len(frequencies) - 1
    fq_str = str(np.round(frequencies[n1])) + '_' + str(np.round(frequencies[n22]))
    datestr = st[0].stats.starttime.strftime('%y%m%d%H')

    if config.plot_results == 'True' or\
            (config.plot_results == 'trigger_only' and trigger is not None):
        bp_plot(config, stack_grid,
                coord_eq, t_begin, t_end, datestr, fq_str,
                coord_sta, st, st_CF,
                frequencies, n1, n22, trigger, arrival_times, Mtau)

    return trigger
