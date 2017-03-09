# -*- coding: utf8 -*-
import numpy as np
from obspy import Stream
from mod_filter_picker import MBfilter_CF, GaussConv


def summary_cf(config, st, rec_memory=None):

    # Compute decay constants
    if config.rosenberger_decay_const is not None:
        rosenberger_decay_const = config.rosenberger_decay_const
    else:
        rosenberger_decay_const = config.decay_const
    # TODO: so far fixed can be added to control file as a separate variable
    sigma_gauss = int(config.decay_const/config.delta/2)

    st_CF = Stream()
    for station in config.stations:
        for wave_type in config.wave_type:
            st_select = st.select(station=station)
            tr_CF = st.select(station=station)[0].copy()
            tr_CF.stats.channel = wave_type
            st_CF.append(tr_CF)
            hos_sigma = config['hos_sigma_' + wave_type]
            MBfilter_CF_kwargs = {
                'st': st_select,
                'frequencies': config.frequencies,
                'CN_HP': config.CN_HP,
                'CN_LP': config.CN_LP,
                'filter_norm': config.filter_norm,
                'filter_npoles': config.filter_npoles,
                'var_w': config.win_type,
                'CF_type': config.ch_function,
                'CF_decay_win': config.decay_const,
                'hos_order': config.hos_order,
                'rosenberger_decay_win': rosenberger_decay_const,
                'rosenberger_filter_power': config.rosenberger_filter_power,
                'rosenberger_filter_threshold':
                    config.rosenberger_filter_threshold,
                'rosenberger_normalize_each':
                    config.rosenberger_normalize_each,
                'wave_type': wave_type,
                'hos_sigma': hos_sigma[station],
                'rec_memory': rec_memory
            }
            HP2, CF, Tn2, Nb2 = MBfilter_CF(**MBfilter_CF_kwargs)

            if config.ch_function == 'envelope':
                tr_CF.data = np.sqrt(np.power(CF, 2).mean(axis=0))
            if config.ch_function == 'kurtosis':
                kurt_argmax = np.amax(CF, axis=0)
                tr_CF.data = GaussConv(kurt_argmax, sigma_gauss)

    #-----resampling CF if wanted-------------------------------------------
    fs_data = st[0].stats.sampling_rate
    if config.sampl_rate_cf:
        if config.sampl_rate_cf < fs_data:
            st_CF.resample(config.sampl_rate_cf)
    else:
        config.sampl_rate_cf = fs_data

    return st_CF


def empty_cf(config, st):
    # creates an empty CF stream between starttime and start_t
    st_CF = Stream()
    for station in config.stations:
        for wave_type in config.wave_type:
            tr_CF = st.select(station=station).copy()
            tr_CF = tr_CF.trim(config.starttime,
                               config.starttime + config.start_t)
            tr_CF = tr_CF[0]
            tr_CF.data = np.zeros(tr_CF.data.shape)
            tr_CF.stats.channel = wave_type
            st_CF.append(tr_CF)

    #-----resampling CF if wanted-------------------------------------------
    fs_data = st[0].stats.sampling_rate
    if config.sampl_rate_cf:
        if config.sampl_rate_cf < fs_data:
            # we don't need to resample if there are less than 2 points
            if len(st_CF[0]) < 2:
                for tr_CF in st_CF:
                    tr_CF.stats.sampling_rate = config.sampl_rate_cf
            else:
                st_CF.resample(config.sampl_rate_cf)
    else:
        config.sampl_rate_cf = fs_data

    return st_CF
