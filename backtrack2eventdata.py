#!/usr/bin/env python
# -*- coding: utf8 -*-
import sys
import os
from tatka_modules.bp_types import Trigger, Pick
from obspy import read
from obspy.core import AttribDict

"""
using output file xxx_OUT2_grouped.dat of BackProj.py
to create events data
to run: ./backtrack2eventdata.py xxx_OUT2_grouped.dat (station.dat --> optional)
"""

def main():
    main_folder = 'Eventdata'
    basepath_data = 'example_data'
    out_data_format = "SAC"

    # data length to cut before and after the P-wave arrival (in seconds)
    pre_P = 20.0
    post_P = 150.0
    ##
    if len(sys.argv) == 1:
        print "this_code  <triggered_location xxx_OUT2_grouped.dat>"
        sys.exit(1)
    else:
        infile = sys.argv[1]
    ##
    if len(sys.argv) > 2:
        station_file = sys.argv[2]
    if len(sys.argv) == 2:
        station_file = None

    ##--reading station information
    coord_sta = {}
    if station_file:
        for line in open(station_file, 'r'):
            data = line.split()
            coord_sta[data[0]] = map(float, data[2:5])

    ## -----reading output file saving parameters as trigger type--------------------------------------
    triggers = []
    for line in open(infile, 'r'):
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

    ##---creating output directories if they do not exist----------------------------------------------
    if not os.path.exists(main_folder):
        os.mkdir(main_folder)
    for trigger in triggers:
        out_event_folder = os.path.join(main_folder, trigger.eventid)
        if not os.path.exists(out_event_folder):
            os.mkdir(out_event_folder)

    ##---outputting event files and data--------------------------------------------------------------
    for trigger in triggers:

    ##-------writting enevtid.dat file ---------------------------------------------------------------
        event_dat_base = '.'.join((trigger.eventid,'dat'))
        event_dat = os.path.join(out_event_folder,event_dat_base)
        with open(event_dat,'w') as f:
            f.write(str(trigger) + '\n')
            for pick in trigger.picks:
                f.write(str(pick) + '\n')

    ##-------writting enevtid_nll.dat file -----------------------------------------------------------
        event_dat_base2 = '.'.join((trigger.eventid+'_nll','dat'))
        event_dat2 = os.path.join(out_event_folder,event_dat_base2)
        with open(event_dat2,'w') as f:
            f.write('#%s %f %f %f %s\n' % (trigger.eventid, trigger.lon, trigger.lat, trigger.z, trigger.origin_time))
            for pick in trigger.picks:
                f.write('%-6s ?    ?    ? %-6s ? ' % (pick.station, pick.arrival_type))
                time = trigger.origin_time + pick.pick_time
                f.write('%s ' % time.strftime('%Y%m%d'))
                f.write('%s ' % time.strftime('%H%M'))
                f.write('%s.' % time.strftime('%S'))
                msec = int(round(int(time.strftime('%f'))/100.))
                f.write('%04d ' % msec)
                f.write('GAU  5.00e-02  0.00e+00  0.00e+00  0.00e+00')
                f.write('\n')

    ## reading data, cutting events and saving data in specified format-------------------------------
        for pick in trigger.picks:
            if pick.arrival_type is not 'P':
                continue

            read_starttime = trigger.origin_time + pick.theor_time - pre_P
            read_endtime = trigger.origin_time + pick.theor_time + post_P
            filename = os.path.join(basepath_data, '*' + pick.station + '*')

            st = read(filename, starttime=read_starttime,
                      endtime=read_endtime)

            for tr in st:
                file_out_base = '.'.join((trigger.eventid,
                                         tr.stats.station,
                                         tr.stats.channel,
                                         tr.stats.network,
                                         tr.stats.location,
                                          out_data_format))
                file_out_data = os.path.join(out_event_folder, file_out_base)
                if out_data_format == 'SAC':
                    tr.stats.sac = AttribDict()
                    tr.stats.sac.kevnm = str(trigger.eventid)
                    if station_file:
                        tr.stats.sac.stla = coord_sta[tr.stats.station][0]
                        tr.stats.sac.stlo = coord_sta[tr.stats.station][1]
                        tr.stats.sac.stel = coord_sta[tr.stats.station][2]

                tr.write(file_out_data, format=out_data_format)


if __name__ == '__main__':
    main()
