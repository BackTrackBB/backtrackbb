# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import scipy as sp
from .LocalCC import LocalCC


def sta_GRD_Proj(args):
    return _sta_GRD_Proj(*args)


def _sta_GRD_Proj(config, sta_wave1, sta_wave2, sig1, sig2, t_b, tau_max=None):

    max_lag = config.time_lag
    if tau_max is not None:
        max_lag = tau_max

    start_time = config.starttime

    if config.do_smooth_lcc:
        sigma = config.smooth_lcc
    else:
        sigma = None

    if config.sampl_rate_cf:
        samp_rate = config.sampl_rate_cf
    else:
        samp_rate = config.sampl_rate

    t_lag, local_cc, arrival1, arrival2 =\
        LocalCC(sig1, sig2, samp_rate, max_lag, start_time+t_b, sigma)

    # Max value of local_cc in given window
    local_cc_1d = np.amax(local_cc, axis=1)

    # Projecting LCC for the station pair on the grid of theoretical t_times
    # Check which of the functions for the interpolations is faster?
    function = sp.interpolate.UnivariateSpline(t_lag, local_cc_1d, k=1, s=0)

    return function, arrival1, arrival2, sta_wave1, sta_wave2
