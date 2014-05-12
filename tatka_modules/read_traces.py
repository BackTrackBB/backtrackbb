# -*- coding: utf8 -*-
import os
from glob import glob
from obspy.core import read, Stream

def read_traces(config):
    basepath = config.data_dir
    if config.data_day:
        basepath = os.path.join(basepath, config.data_day)
        if config.data_hours:
            basepath = os.path.join(basepath, config.data_hours)

    kwargs = {}
    if config.data_format:
        kwargs['format'] = config.data_format

    tmpst = Stream()
    for filename in glob(os.path.join(basepath, '*')):
        try:
            tmpst += read(filename, **kwargs)
        except Exception:
            continue

    # Get the intersection between the list of available stations
    # and the list of requested stations:
    tmpst_stations = [tr.stats.station for tr in tmpst.select(channel=config.channel)]
    stations = sorted(set(tmpst_stations) & set(config.stations))

    # Retain only requested channel and stations:
    st = Stream(tr for tr in tmpst.select(channel=config.channel)
                if tr.stats.station in stations)

    print 'Number of traces in stream = ', len(st)

    st.sort()
    return st, stations
