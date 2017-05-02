# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import os
import numpy as np
import itertools
from collections import defaultdict
from multiprocessing import Pool
from scipy.ndimage.interpolation import zoom
from .bp_types import Trigger
from .map_project import rect2latlon
from .grid_projection import sta_GRD_Proj
from .plot import bp_plot
from .mod_bp_TrigOrig_time import TrOrig_time
from .NLLGrid import NLLGrid
from .summary_cf import summary_cf


def init_worker():
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def slice_indexes(i, j, k, si, sj, sk, ni, nj, nk):
    i1 = i - si if i - si > 0 else 0
    i2 = i + si if i + si < ni else ni
    j1 = j - sj if j - sj > 0 else 0
    j2 = j + sj if j + sj < nj else nj
    k1 = k - sk if k - sk > 0 else 0
    k2 = k + sk if k + sk < nk else nk
    return i1, i2, j1, j2, k1, k2


def run_btbb(args):
    return _run_btbb(*args)


def _run_btbb(config, st, st_CF, t_begin,
              coord_sta, GRD_sta, coord_eq,
              rec_memory=None, async_plotter=None):

    t_end = t_begin + config.time_lag

    # Create stack_grid using a time grid as model
    gr = list(list(GRD_sta.values())[0].values())[0]
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

    if rec_memory is not None:
        st_cut = st.copy()
        st_cut = st_cut.trim(config.starttime + t_begin,
                             config.starttime + t_end)
        st_CF_cut = summary_cf(config, st_cut, rec_memory=rec_memory)
        st_CF += st_CF_cut
        st_CF.merge(method=1)
    else:
        st_CF_cut = st_CF.copy()
        st_CF_cut = st_CF_cut.trim(config.starttime + t_begin,
                                   config.starttime + t_end)

    if config.ignore_noisy_CF:
        sums = [('%s' % str(tr.stats.station), '%s' % str(tr.stats.channel),
                (tr.data/tr.max()).sum()) for tr in st_CF_cut]

        min_sum = min(s[1] for s in sums)

        # A noisy CF has an integral at least two times larger than
        # the smaller one
        # TODO: parametrize?
        noisy_sta_wave = [(s[0], s[1]) for s in sums if s[2] >= 2*min_sum]
    else:
        noisy_sta_wave = []

    # Prepare the arglist for sta_Grd_Proj()
    Mtau = []
    arrival_times = defaultdict(dict)
    arglist = []
    sta_wave = [(sta, wave)
                for wave in config.wave_type for sta in config.stations]

    # Remove noisy CFs, if this leaves us with at least 4 CFs (3 sta are few)
    # (otherwhise use all the CFs)
    if (len(sta_wave) - len(noisy_sta_wave)) >= 4:
        sta_wave = [s for s in sta_wave if s not in noisy_sta_wave]
    for sta_wave1, sta_wave2 in itertools.combinations(sta_wave, 2):
        sta1 = sta_wave1[0]
        sta2 = sta_wave2[0]
        wave1 = sta_wave1[1]
        wave2 = sta_wave2[1]
        arrival_times[sta1][wave1] = []
        arrival_times[sta2][wave2] = []

        x_sta1, y_sta1 = coord_sta[sta1]
        x_sta2, y_sta2 = coord_sta[sta2]
        distance = np.sqrt((x_sta1-x_sta2)**2 + (y_sta1-y_sta2)**2)
        if distance > config.maxSTA_distance:
            continue

        if config.varWin_stationPair:
            tau_max = GRD_sta[sta1][wave1].get_value(
                GRD_sta[sta2][wave1].sta_x,
                GRD_sta[sta2][wave1].sta_y,
                GRD_sta[sta2][wave1].sta_z)
            Mtau.append(np.round(tau_max, 1))
            t_end = t_begin + np.round(tau_max, 1)
        else:
            Mtau = None
            tau_max = None

        trace1 = st_CF_cut.select(station=sta1, channel=wave1)[0]
        trace2 = st_CF_cut.select(station=sta2, channel=wave2)[0]
        sig1 = trace1.data/max(abs(trace1.data))
        sig2 = trace2.data/max(abs(trace2.data))
        len1 = sig1.size
        len2 = sig2.size
        if 0 < abs(len1-len2) < 3:
            min_len = min(len1, len2)
            sig1 = sig1[:min_len]
            sig2 = sig2[:min_len]

        arglist.append((config, sta_wave1, sta_wave2,
                        sig1, sig2, t_begin, tau_max))

    # Run sta_Grd_Proj(), possibly in parallel
    if config.recursive_memory and config.ncpu > 1:
        # When using recursive_memory, parallelization is
        # done here
        rm_pool = Pool(config.ncpu, init_worker)
        try:
            # we need to use map_async() (with a very long timeout)
            # due to a python bug
            # (http://stackoverflow.com/questions/1408356/
            #  keyboard-interrupts-with-pythons-multiprocessing-pool)
            outputs = rm_pool.map_async(sta_GRD_Proj, arglist).get(9999999)
        except KeyboardInterrupt:
            rm_pool.terminate()
            rm_pool.join()
            print('')
            print('Aborting.')
            sys.exit()
        rm_pool.close()
        rm_pool.join()
    else:
        outputs = list(map(sta_GRD_Proj, arglist))

    # Parse outputs and update stack_grid
    for out in outputs:
        proj_function, arrival1, arrival2, sta_wave1, sta_wave2 = out
        sta1 = sta_wave1[0]
        sta2 = sta_wave2[0]
        wave1 = sta_wave1[1]
        wave2 = sta_wave2[1]
        arrival_times[sta1][wave1].append(arrival1)
        arrival_times[sta2][wave2].append(arrival2)
        delta_tt_array =\
            GRD_sta[sta2][wave2].array - GRD_sta[sta1][wave1].array
        stack_grid.array += proj_function(
            delta_tt_array.flatten()).reshape(delta_tt_array.shape)
    stack_grid.array /= len(outputs)
    stack_grid.array **= config.grid_power

    i_max, j_max, k_max = stack_grid.get_ijk_max()

    do_trigger = False
    if config.trigger is not None:
        if stack_grid.max() >= config.trigger:
            do_trigger = True
            trigger_level = config.trigger
    if config.trigger_ellipsoid_max_axis is not None and not do_trigger:
        gp_ell = config.grid_power_ellipsoid
        gp = config.grid_power
        if gp_ell != gp:
            pwr = gp_ell/gp
            stack_grid_pwr = stack_grid.copy()
            stack_grid_pwr.array **= pwr
            ell = stack_grid_pwr.get_xyz_ellipsoid()
            stack_grid.ellipsoid = ell
            stack_grid.xyz_mean = stack_grid_pwr.get_xyz_mean()
        else:
            ell = stack_grid.get_xyz_ellipsoid()
        max_ax = config.trigger_ellipsoid_max_axis
        if (ell.len1 <= max_ax and
                stack_grid.max() >= config.trigger_ellipsoid):
            do_trigger = True
            trigger_level = config.trigger_ellipsoid
    if config.trigger_probability_range is not None and not do_trigger:
        rng = config.trigger_probability_range
        i_range = int(rng / stack_grid.dx)
        j_range = int(rng / stack_grid.dy)
        k_range = int(rng / stack_grid.dz)
        mask = np.zeros(stack_grid.array.shape)
        ni, nj, nk = mask.shape
        idx = slice_indexes(i_max, j_max, k_max,
                            i_range, j_range, k_range,
                            ni, nj, nk)
        stack_grid.box_idx = idx
        i1, i2, j1, j2, k1, k2 = idx
        mask[i1:i2, j1:j2, k1:k2] = 1
        mask_array = stack_grid.array * mask
        prob = mask_array.sum() / stack_grid.array.sum()
        if prob >= config.trigger_probability:
            do_trigger = True
            trigger_level =\
                0.5 * (stack_grid.array.max() - stack_grid.array.min())

    if config.cut_data:
        start_tw = config.cut_start + t_begin
    else:
        start_tw = t_begin
        config.cut_start = 0.

    trigger = None
    if do_trigger:
        if config.max_subdivide is not None:
            # We "zoom" the grid in a region close to the maximum
            zoom_factor = float(config.max_subdivide)
            slice_factor = 4  # this is hardcoded, for now
            sf = int(slice_factor)
            ni, nj, nk = stack_grid.array.shape
            idx = slice_indexes(i_max, j_max, k_max,
                                sf, sf, sf,
                                ni, nj, nk)
            i1, i2, j1, j2, k1, k2 = idx
            slice_grid = stack_grid[i1:i2, j1:j2, k1:k2]
            zoom_slice_grid = zoom(slice_grid, zoom_factor)
            # check if zoom_slice_grid is not empty.
            # TODO: why can it be empty?
            if zoom_slice_grid.size > 0:
                zoom_i_max, zoom_j_max, zoom_k_max =\
                    [a[0]/zoom_factor
                     for a in
                     np.where(zoom_slice_grid == np.max(zoom_slice_grid))]
                i_max = zoom_i_max + i_max - sf
                i_max = i_max if i_max > 0 else 0
                j_max = zoom_j_max + j_max - sf
                j_max = j_max if j_max > 0 else 0
                k_max = zoom_k_max + k_max - sf
                k_max = k_max if k_max > 0 else 0

        trigger = Trigger()
        trigger.trigger_level = trigger_level
        trigger.x, trigger.y, trigger.z =\
            stack_grid.get_xyz(i_max, j_max, k_max)
        trigger.i, trigger.j, trigger.k = i_max, j_max, k_max
        trigger.max_grid = np.round(stack_grid.max(), 4)
        trigger.ntraces = len(sta_wave)
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
        print(trigger)

    if config.save_projGRID is True or\
            (config.save_projGRID == 'trigger_only' and trigger is not None):
        print('Saving projection grid to file.')
        basename = os.path.join(config.out_dir, 'out_t%05.1f' % t_begin)
        stack_grid.write_hdr_file(basename)
        stack_grid.write_buf_file(basename)

    # Plotting------------------------------------------------------------------
    if config.plot_results == 'True' or\
            (config.plot_results == 'trigger_only' and trigger is not None):

        bp_plot(config, stack_grid,
                coord_eq, t_begin, t_end,
                coord_sta, st, st_CF,
                trigger,
                arrival_times,
                noisy_sta_wave,
                Mtau, async_plotter)

    return trigger
