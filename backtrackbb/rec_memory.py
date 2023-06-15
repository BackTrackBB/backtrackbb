# -*- coding: utf8 -*-
import itertools
from .bp_types import RecursiveMemory


def init_recursive_memory(config):
    n_bands = config.n_freq_bands
    nsamples = int(config.time_lag / config.delta)
    overlap = int(config.t_overlap / config.delta)
    # Create a dictionary of memory objects
    rec_memory = dict()
    for trid, wave in itertools.product(config.trids, config.wave_type):
        # Each entry of the dictionary is a list of memory objects
        # (with n_bands elements)
        rec_memory[(trid, wave)] =\
            [RecursiveMemory(trid=trid, wave=wave, band=n,
                             nsamples=nsamples, overlap=overlap,
                             filter_npoles=config.filter_npoles)
             for n in range(n_bands)]
    return rec_memory
