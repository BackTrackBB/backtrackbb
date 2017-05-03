# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import numpy as np
from multiprocessing import Pool
from ..mod_setup import configure
from ..read_traces import read_traces
from ..init_filter import init_filter
from ..read_grids import read_grids
from ..summary_cf import summary_cf, empty_cf
from ..map_project import get_transform
from ..mod_utils import read_locationTremor, read_locationEQ
from ..plot import plt_SummaryOut
from ..rec_memory import init_recursive_memory
from ..mod_btbb import run_btbb
from ..AsyncPlotter import AsyncPlotter

DEBUG = False


def main():
    config = configure()

    var_twin = config.varWin_stationPair
    print('use of var time window for location:', var_twin)
    #---Reading data---------------------------------------------------------
    st = read_traces(config)
    #------------------------------------------------------------------------
    t_bb = np.arange(config.start_t, config.end_t,
                     config.time_lag - config.t_overlap)
    # selecting the time windows that do not exceed the length of the data---
    t_ee = t_bb + config.time_lag
    data_length = st[0].stats.endtime - st[0].stats.starttime
    t_bb = t_ee[t_ee <= data_length] - config.time_lag
    print('Number of time windows = ', len(t_bb))

    loc_infile = None
    location_jma = None
    if config.catalog_dir:
        if config.data_day:
            loc_infile = os.path.join(config.catalog_dir,
                                      config.data_day+config.tremor_file)
        location_jma = os.path.join(config.catalog_dir, config.eq_file)
    #------------------------------------------------------------------------

    #--Read grids of theoretical travel-times--------------------------------
    GRD_sta, coord_sta = read_grids(config)

    #---remove mean and trend------------------------------------------------
    st.detrend(type='constant')
    st.detrend(type='linear')

    #---init filtering parameters
    init_filter(config)

    if config.recursive_memory:
        rec_memory = init_recursive_memory(config)
        st_CF = empty_cf(config, st)
        if DEBUG:
            st_CF2 = summary_cf(config, st)
    else:
        rec_memory = None
        st_CF = summary_cf(config, st)

    #---Take the first grid as reference ------------------------------------
    grid1 = list(list(GRD_sta.values())[0].values())[0]

    #---init map projection--------------------------------------------------
    if grid1.proj_name:
        get_transform(grid1.proj_name,
                      grid1.orig_lat, grid1.orig_lon,
                      grid1.first_std_paral,
                      grid1.second_std_paral,
                      grid1.map_rot,
                      grid1.proj_ellipsoid)

    #----geographical coordinates of the eq's epicenter----------------------
    coord_eq = None
    if loc_infile:
        coord_eq = read_locationTremor(loc_infile, config.data_hours,
                                       config.lat_orig, config.lon_orig)
    coord_jma = None
    if location_jma:
        coord_jma = read_locationEQ(location_jma,
                                    config.data_day, config.data_hours,
                                    config.lat_orig, config.lon_orig)
    #------------------------------------------------------------------------
    print('starting BPmodule')

    # Create out_dir, if it doesn't exist
    if not os.path.exists(config.out_dir):
        os.mkdir(config.out_dir)

    datestr = st[0].stats.starttime.strftime('%y%m%d%H')
    fq_str = '%s_%s' % (np.round(config.frequencies[0]),
                        np.round(config.frequencies[-1]))
    ch_str = str(config.channel)[1:-1].replace("'", "")
    file_out_base = '_'.join((
        datestr,
        str(len(config.frequencies)) + 'fq' + fq_str + 'hz',
        str(config.decay_const) + str(config.sampl_rate_cf) +
        str(config.smooth_lcc) + str(config.t_overlap),
        config.ch_function,
        ch_str,
        ''.join(config.wave_type),
        'trig' + str(config.trigger)
        ))

    file_out_triggers = file_out_base + '_OUT2.dat'
    file_out_triggers = os.path.join(config.out_dir, file_out_triggers)

    file_out_fig = file_out_base + '_FIG2.' + config.plot_format
    file_out_fig = os.path.join(config.out_dir, file_out_fig)

    #---running program------------------------------------------------------
    if config.ncpu > 1 and config.recursive_memory \
            and config.plot_results != 'False':
        global async_plotter
        async_plotter = AsyncPlotter()
    else:
        async_plotter = None

    arglist = [
               (config,
                st, st_CF, t_begin,
                coord_sta, GRD_sta, coord_eq,
                rec_memory, async_plotter)
               for t_begin in t_bb
              ]
    print('Running on %d thread%s' % (config.ncpu, 's' * (config.ncpu > 1)))
    if config.ncpu > 1 and not config.recursive_memory:
        # parallel execution
        global pool
        pool = Pool(config.ncpu, init_worker)
        # we need to use map_async() (with a very long timeout)
        # due to a python bug
        # (http://stackoverflow.com/questions/1408356/
        #  keyboard-interrupts-with-pythons-multiprocessing-pool)
        p_outputs = pool.map_async(run_btbb, arglist).get(9999999)
        pool.close()
        pool.join()
    else:
        # serial execution
        # (but there might be a parallelization step
        # inside run_btbb, if we're using recursive_memory)
        p_outputs = map(run_btbb, arglist)
    triggers = list(filter(None, p_outputs))

    if async_plotter is not None:
        async_plotter.join()

    #----------Outputs-------------------------------------------------------
    #writing output
    eventids = []
    with open(file_out_triggers, 'w') as f:
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

    #-plot summary output----------------------------------------------------
    plt_SummaryOut(config, grid1, st_CF, st, coord_sta,
                   triggers, t_bb, datestr,
                   coord_eq, coord_jma, file_out_fig)

    if DEBUG and config.recursive_memory:
        import matplotlib.pyplot as plt
        CF = st_CF.select(station=config.stations[0])[0]
        CF2 = st_CF2.select(station=config.stations[0])[0]
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.plot(CF, linewidth=2)
        ax1.plot(CF2[0:len(CF)])
        ax2 = fig.add_subplot(212, sharex=ax1)
        ax2.plot(CF - CF2[0:len(CF)])
        plt.show()


def init_worker():
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def run():
    try:
        pool = None
        async_plotter = None
        main()
    except KeyboardInterrupt:
        if pool is not None:
            pool.terminate()
            pool.join()
        if async_plotter is not None:
            async_plotter.terminate()
        print('')
        print('Aborting.')


if __name__ == '__main__':
    run()
