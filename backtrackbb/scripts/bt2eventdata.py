# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
from obspy.core import AttribDict
from ..mod_setup import configure
from ..bp_types import Trigger, Pick
from ..read_traces import read_traces

"""
using output file xxx_OUT2_grouped.dat of btbb.py
to create events data
to run: ./bt2eventdata.py config_file xxx_OUT2_grouped.dat
        (station.dat --> optional)
"""


def main():
    config = configure('backtrack2eventdata')

    ##--reading station information
    coord_sta = {}
    if config.options.station_file:
        for line in open(config.options.station_file, 'r'):
            data = line.split()
            coord_sta[data[0]] = map(float, data[2:5])

    # -----reading output file saving parameters as trigger type-----
    triggers = []
    for line in open(config.options.trigger_file, 'r'):
        try:
            trigger = Trigger()
            trigger.from_str(line)
            triggers.append(trigger)
        except ValueError:
            try:
                pick = Pick()
                pick.from_str(line)
                triggers[-1].add_pick(pick)
            except ValueError:
                continue

    # -----creating output directories if they do not exist-----
    if not os.path.exists(config.event_dir):
        os.mkdir(config.event_dir)

    # -----Reading data-----
    st = read_traces(config)

    for trigger in triggers:
        if trigger.ntraces < config.min_ntraces:
            continue
        if trigger.max_grid < config.min_trigger:
            continue
        print(trigger.eventid)
        out_event_dir = os.path.join(config.event_dir, trigger.eventid)
        if not os.path.exists(out_event_dir):
            os.mkdir(out_event_dir)

        # -----writing eventid.dat file-----
        event_dat_base = '.'.join((trigger.eventid, 'dat'))
        event_dat = os.path.join(out_event_dir, event_dat_base)
        with open(event_dat, 'w') as f:
            f.write(str(trigger) + '\n')
            for pick in trigger.picks:
                f.write(str(pick) + '\n')

        # -----writing eventid_nll.obs file-----
        event_dat_base = '.'.join((trigger.eventid, 'pick'))
        event_dat = os.path.join(out_event_dir, event_dat_base)

        with open(event_dat, 'w') as f:
            for pick in trigger.picks:
                f.write('%-6s ?    ?    ? %-6s ? ' %
                        (pick.station, pick.arrival_type))
                if pick.arrival_type == 'P':
                    time = trigger.origin_time + pick.pick_time
                else:
                    time = trigger.origin_time + pick.theor_time
                f.write('%s ' % time.strftime('%Y%m%d'))
                f.write('%s ' % time.strftime('%H%M'))
                f.write('%s.' % time.strftime('%S'))
                msec = int(round(int(time.strftime('%f'))/100.))
                f.write('%04d ' % msec)
                # We approximate the error with decay_const/10
                # TODO: improve this?
                decay_const = float(config.decay_const)
                f.write('GAU  %.2e  0.00e+00  0.00e+00  0.00e+00' %
                        (decay_const/10.))
                f.write('\n')

        # -----cutting events and saving data in specified format-----
        for pick in trigger.picks:
            if pick.arrival_type != 'P':
                continue

            trim_starttime = (trigger.origin_time + pick.theor_time
                              - config.pre_P)
            trim_endtime = (trigger.origin_time + pick.theor_time
                            + config.post_P)

            st_select = st.copy().select(station=pick.station)
            st_select.trim(trim_starttime, trim_endtime)

            for tr in st_select:
                file_out_base = '.'.join((trigger.eventid,
                                         tr.stats.station,
                                         tr.stats.channel,
                                         tr.stats.network,
                                         tr.stats.location,
                                         config.out_data_format))
                file_out_data = os.path.join(out_event_dir, file_out_base)
                if config.out_data_format == 'sac':
                    tr.stats.sac = AttribDict()
                    tr.stats.sac.kevnm = str(trigger.eventid)
                    if config.options.station_file:
                        tr.stats.sac.stla = coord_sta[tr.stats.station][0]
                        tr.stats.sac.stlo = coord_sta[tr.stats.station][1]
                        tr.stats.sac.stel = coord_sta[tr.stats.station][2]

                tr.write(file_out_data, format=config.out_data_format)


if __name__ == '__main__':
    main()
