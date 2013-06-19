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
    
    st = Stream()
    for filename in glob(os.path.join(basepath, '*')):
        try:
            tmpst = read(filename, **kwargs)
        except Exception:
            continue

        # Retain only requested component:
        tmpst = tmpst.select(component=config.component)
        # Retain only requested stations:
        for tr in tmpst:
            if tr.stats.station in config.stations:
                st.append(tr)
                break

    print 'Number of stations in stream = ', len(st)

    st.sort()
    return st
