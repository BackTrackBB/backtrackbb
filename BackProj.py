#!/usr/bin/env python
# -*- coding: utf8 -*-
import sys
import os
import copy
import numpy as np
from tatka_modules.read_traces import read_traces
from tatka_modules.mod_filter_picker import make_LinFq, make_LogFq
from tatka_modules.read_grids import read_grids
from tatka_modules.summary_cf import summary_cf, empty_cf
from tatka_modules.map_project import get_transform
from tatka_modules.mod_utils import read_locationTremor, read_locationEQ
from tatka_modules.mod_groupe_trigs import groupe_triggers
from tatka_modules.plot import plt_SummaryOut
from tatka_modules.parse_config import parse_config
from tatka_modules.rec_memory import init_recursive_memory
from tatka_modules.mod_backproj import run_BackProj
from multiprocessing import Pool

DEBUG = False

def main():
    if len(sys.argv) != 2:
        print "this_code  <input_config_file>"
        sys.exit(1)
    else:
        Config_file = sys.argv[1]
    if not os.path.isfile(Config_file):
        print "File {0} does not exist".format(Config_file)
        sys.exit(1)

    #---Input parameters for BProj run----------------------------------------
    config = parse_config(Config_file)
    var_twin = config.varWin_stationPair
    print 'use of var time window for location:', var_twin
    #---Reading data---------------------------------------------------------
    st = read_traces(config)
    #------------------------------------------------------------------------
    t_bb = np.arange(config.start_t, config.end_t, config.time_lag - config.t_overlap)
    # selecting the time windows that do not exceed the length of the data---
    t_ee = t_bb + config.time_lag
    data_length = st[0].stats.endtime - st[0].stats.starttime
    t_bb = t_ee[t_ee<=data_length] - config.time_lag
    print 'Number of time windows = ', len(t_bb)

    loc_infile = None
    location_jma = None
    if config.catalog_dir:
        if config.data_day:
            loc_infile = os.path.join(config.catalog_dir, config.data_day+config.tremor_file)
        location_jma = os.path.join(config.catalog_dir, config.eq_file)
    #------------------------------------------------------------------------

    #--Read grids of theoretical travel-times--------------------------------
    GRD_sta, coord_sta = read_grids(config)

    #--- cut the data to the selected length dt------------------------------
    if config.cut_data:
        st.trim(st[0].stats.starttime+config.cut_start,
                st[0].stats.starttime+config.cut_start+config.cut_delta)
    else:
        config.cut_start = 0.

    #---remove mean and trend------------------------------------------------
    st.detrend(type='constant')
    st.detrend(type='linear')

    #---Calculating frequencies for MBFilter---------------------------------
    if config.band_spacing == 'lin':
        frequencies = make_LinFq(config.f_min, config.f_max, config.delta, config.n_freq_bands)
    elif config.band_spacing == 'log':
        frequencies = make_LogFq(config.f_min, config.f_max, config.delta, config.n_freq_bands)
    n1 = 0
    n2 = len(frequencies)
    n22 = len(frequencies) - 1
    print 'frequencies for filtering in (Hz):', frequencies[n1:n2]

    if config.recursive_memory:
        rec_memory = init_recursive_memory(config)
        st_CF = empty_cf(config, st)
        if DEBUG:
            st_CF2 = summary_cf(config, st, frequencies)
    else:
        rec_memory = None
        st_CF = summary_cf(config, st, frequencies)

    #---Take the first grid as reference ------------------------------------
    grid1 = GRD_sta.values()[0].values()[0]

    #---init map projection--------------------------------------------------
    if grid1.proj_name!='NONE':
        get_transform(grid1.proj_name,
                      grid1.orig_lat, grid1.orig_lon,
                      grid1.first_std_paral,
                      grid1.second_std_paral,
                      grid1.map_rot,
                      grid1.ellipsoid)

    #----geographical coordinates of the eq's epicenter----------------------
    coord_eq = None
    if loc_infile:
        coord_eq = read_locationTremor(loc_infile, config.data_hours,
                                           config.lat_orig, config.lon_orig)
    coord_jma = None
    if location_jma:
        coord_jma = read_locationEQ(location_jma, config.data_day, config.data_hours,
                                        config.lat_orig, config.lon_orig)
    #------------------------------------------------------------------------
    print 'starting BPmodule'

    # Create out_dir, if it doesn't exist
    if not os.path.exists(config.out_dir):
        os.mkdir(config.out_dir)

    datestr = st[0].stats.starttime.strftime('%y%m%d%H')
    fq_str = str(np.round(frequencies[n1])) + '_' + str(np.round(frequencies[n22]))
    file_out_base = '_'.join((
        datestr,
        str(len(frequencies)) + 'fq' + fq_str + 'hz',
        str(config.decay_const) + str(config.sampl_rate_cf) + str(config.smooth_lcc) + str(config.t_overlap),
        config.ch_function,
        config.channel,
        ''.join(config.wave_type),
        'trig' + str(config.trigger)
        ))

    file_out_data = file_out_base + '_OUT2.dat'
    file_out_data = os.path.join(config.out_dir, file_out_data)

    file_out_data2 = file_out_base + '_OUT2_grouped.dat'
    file_out_data2 = os.path.join(config.out_dir, file_out_data2)

    file_out_fig = file_out_base + '_FIG2.' + config.plot_format
    file_out_fig = os.path.join(config.out_dir, file_out_fig)

    #---running program------------------------------------------------------
    if config.recursive_memory:
        # We use plot_pool to run asynchronus plotting in run_Backproj
        plot_pool = Pool(config.ncpu)
    else:
        plot_pool = None

    arglist = [
               (config,
                st, st_CF, t_begin, frequencies,
                coord_sta, GRD_sta,
                coord_eq, rec_memory, plot_pool)
               for t_begin in t_bb
              ]
    print 'Running on %d thread%s' % (config.ncpu, 's' * (config.ncpu > 1))
    if config.ncpu > 1 and not config.recursive_memory:
        # parallel execution
        p = Pool(config.ncpu) #defining number of jobs
        p_outputs = p.map(run_BackProj, arglist)
        p.close() #no more tasks
        p.join() #wrap up current tasks
    else:
        # serial execution
        # (but there might be a parallelization step
        # inside run_BackProj, if we're using recursive_memory)
        p_outputs = map(run_BackProj, arglist)
    triggers = filter(None, p_outputs)

    if config.recursive_memory:
        plot_pool.close()

    #----------Outputs-------------------------------------------------------
    #writing output
    eventids = []
    with open(file_out_data,'w') as f:
        for trigger in triggers:
            # check if eventid already exists
            while trigger.eventid in eventids:
                # increment the last letter by one
                evid = list(trigger.eventid)
                evid[-1] = chr(ord(evid[-1]) + 1)
                trigger.eventid = ''.join(evid)
            eventids.append(trigger.eventid)
            f.write(str(trigger) + '\n')
            # sort picks by station
            picks = sorted(trigger.picks, key=lambda x: x.station)
            for pick in picks:
                f.write(str(pick) + '\n')
    if triggers:
        #----sorting triggers----and grouping triggered locations----------------
        sorted_trigs = copy.copy(triggers)
        sorted_trigs = groupe_triggers(sorted_trigs,config)

        #writing sorted triggers in a second output file-------------------------
        #writing output
        eventids = []
        with open(file_out_data2,'w') as f:
            for trigger in sorted_trigs:
                # check if eventid already exists
                while trigger.eventid in eventids:
                    # increment the last letter by one
                    evid = list(trigger.eventid)
                    evid[-1] = chr(ord(evid[-1]) + 1)
                    trigger.eventid = ''.join(evid)
                eventids.append(trigger.eventid)
                f.write(str(trigger) + '\n')
                # sort picks by station
                picks = sorted(trigger.picks, key=lambda x: x.station)
                for pick in picks:
                    f.write(str(pick) + '\n')

        #-plotting output-------------------------------------------------------
        plt_SummaryOut(config, grid1, st_CF, st, coord_sta,
                       sorted_trigs, t_bb, datestr, frequencies[n1], frequencies[n22],
                       coord_eq, coord_jma, file_out_fig)
    else:
    #-plotting output--------------------------------------------------------
        plt_SummaryOut(config, grid1, st_CF, st, coord_sta,
                       triggers, t_bb, datestr, frequencies[n1], frequencies[n22],
                       coord_eq, coord_jma, file_out_fig)

    if DEBUG:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.plot(st_CF[0], linewidth=2)
        ax1.plot(st_CF2[0][0:len(st_CF[0])])
        ax2 = fig.add_subplot(212, sharex=ax1)
        ax2.plot(st_CF[0] - st_CF2[0][0:len(st_CF[0])])
        plt.show()


if __name__ == '__main__':
    main()
