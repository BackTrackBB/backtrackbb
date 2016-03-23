# -*- coding: utf8 -*-
import sys
import os
from glob import glob
from obspy.core import read, Stream,UTCDateTime

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
            if config.start_time:
                tmpst += read(filename,
                              starttime = UTCDateTime(config.start_time),
                              endtime = UTCDateTime(config.end_time),**kwargs)
            else:
                tmpst += read(filename, **kwargs)
        except Exception:
            continue

    # Get the intersection between the list of available stations
    # and the list of requested stations:
    tmpst_select = Stream()
    for ch in config.channel:
        tmpst_select += tmpst.select(channel=ch)
    tmpst_stations = [tr.stats.station for tr in tmpst_select]
    stations = sorted(set(tmpst_stations) & set(config.stations))

    # Retain only requested channel and stations:
    st = Stream(tr for tr in tmpst_select if tr.stats.station in stations)
    if not st:
        print 'Could not read any trace!'
        sys.exit(1)
    st.sort()

    # Check sampling rate
    config.delta = None
    for tr in st:
        tr.detrend(type='constant')
        tr.taper(type='hann', max_percentage=0.005, side='left')
        sampling_rate = tr.stats.sampling_rate
        # Resample data, if requested
        if config.sampl_rate_data:
            if sampling_rate >= config.sampl_rate_data:
                dec_ct = int(sampling_rate/config.sampl_rate_data)
                tr.decimate(dec_ct, strict_length=False, no_filter=True)
            else:
                raise ValueError, 'Sampling frequency for trace %s is lower than %s' % (tr.id, config.sampl_rate_data)
        delta = tr.stats.delta
        if config.delta is None:
            config.delta = delta
        else:
            if delta != config.delta:
                raise ValueError, 'Trace %s has different delta: %s (expected: %s)' % (tr.id, delta, config.delta)
    # Recompute sampling rate after resampling
    config.sampl_rate_data = st[0].stats.sampling_rate

    print 'Number of traces in stream = ', len(st)

    # Check for common starttime and endtime of the traces
    st_starttime = max([tr.stats.starttime for tr in st])
    st_endtime = min([tr.stats.endtime for tr in st])
    if config.start_time:
        st.trim(max(st_starttime, UTCDateTime(config.start_time)),
                min(st_endtime, UTCDateTime(config.end_time)))
    else:
        st.trim(st_starttime, st_endtime)

    #--- cut the data to the selected length dt------------------------------
    if config.cut_data:
        st.trim(st[0].stats.starttime+config.cut_start,
                st[0].stats.starttime+config.cut_start+config.cut_delta)
    else:
        config.cut_start = 0.

    config.starttime = st[0].stats.starttime

    # attach station list and trace ids to config file
    config.stations = stations
    config.trids = [tr.id for tr in st]
    return st
