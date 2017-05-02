from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import scipy as sp
from .rec_cc import local_CCr


def LocalCC(in_signal1, in_signal2, samp_rate, max_time_lag, zero_time,
            sigma=None):

    n_lags = 2 * int(max_time_lag * samp_rate)
    cc_time_lags = sp.linspace(-n_lags/2/samp_rate,
                               n_lags/2/samp_rate,
                               n_lags, endpoint=True)
    cc = local_CCr(in_signal1, in_signal2,
                   max_time_lag, samp_rate, sigma)

    lag, time = np.unravel_index(cc.argmax(), cc.shape)

    time = zero_time + time/samp_rate
    lag = lag/samp_rate - max_time_lag

    arrival1 = time - lag/2.
    arrival2 = time + lag/2.

    return cc_time_lags, cc, arrival1, arrival2
