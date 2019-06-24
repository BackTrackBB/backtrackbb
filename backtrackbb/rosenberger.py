# -*- coding: utf8 -*-
# rosenberger.py
#
# Python interface to C code for P and non-P wavefield separation
# by Andreas Rosenberger.
#
# Rosenberger, A., 2010. Real-Time Ground-Motion Analysis:
# Distinguishing P and S Arrivals in a Noisy Environment.
# Bull. Seismol. Soc. Am. 100, 1252–1262. doi: 10.1785/0120090265
#
# (c) 2014 - Claudio Satriano <satriano@ipgp.fr>,
#            Pierre Romanet <romanet@ipgp.fr>
# (c) 2014 - 2018 Claudio Satriano <satriano@ipgp.fr>
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from ctypes import CDLL, c_int, c_float, c_char, c_void_p
from numpy.ctypeslib import ndpointer
from future.utils import PY2
try:
    from .lib_names import get_lib_path
except (ImportError, ValueError):
    from lib_names import get_lib_path


lib_rosenberger = CDLL(get_lib_path('lib_rosenberger'))
lib_rosenberger.rosenberger.argtypes = [
        ndpointer(dtype=np.float64),  # dataX
        ndpointer(dtype=np.float64),  # dataY
        ndpointer(dtype=np.float64),  # dataZ
        ndpointer(dtype=np.float64),  # pol_filter
        c_int,  # npts
        c_float,  # lambda
        c_float,  # delta
        c_char,  # proj
        c_char  # rl_filter
        ]
lib_rosenberger.rosenberger.restype = c_void_p


def rosenberger(dataX, dataY, dataZ,
                lambda_, delta=1, proj=False, rl_filter=False,
                pol_filter_power=1., pol_filter_threshold=None,
                normalize_each=False):
    """
    Separate P and non-P wavefield from 3-component data.

    Returns two set of 3-component traces and the polarization filter.
    """
    dataX = np.array(dataX, dtype=np.float64)
    dataY = np.array(dataY, dtype=np.float64)
    dataZ = np.array(dataZ, dtype=np.float64)

    # Normalize and multiply by a large factor to avoid small numbers in SVD
    normX = np.abs(dataX).max()
    normY = np.abs(dataY).max()
    normZ = np.abs(dataZ).max()
    if not normalize_each:
        normX = normY = normZ = max(normX, normY, normZ)
    factor = 10000.
    dataX /= normX / factor
    dataY /= normY / factor
    dataZ /= normZ / factor

    pol_filter = np.zeros_like(dataX)

    delta = c_float(delta)
    if PY2:
        proj = chr(proj)
        rl_filter = chr(rl_filter)
    else:
        proj = c_char(proj)
        rl_filter = c_char(rl_filter)

    lib_rosenberger.rosenberger(dataX, dataY, dataZ,
                                pol_filter,
                                len(dataX),
                                lambda_, delta, proj, rl_filter)
    dataX *= normX / factor
    dataY *= normY / factor
    dataZ *= normZ / factor

    pol_filter **= pol_filter_power
    if pol_filter_threshold is not None:
        pol_filter = (pol_filter >= pol_filter_threshold).astype(int)
    data_P = np.vstack((dataZ, dataX, dataY)) * pol_filter[None, :]
    data_S = np.vstack((dataZ, dataX, dataY)) * (1 - pol_filter[None, :])
    return data_P, data_S, pol_filter


def main():
    import math
    from obspy import read
    import matplotlib.pyplot as plt

    # We use the default ObsPy example
    st = read()
    st.filter(type='highpass', freq=1.0)
    maxval = max(np.abs(st.max()))

    time = np.arange(len(st[1].data)) * st[1].stats.delta

    # Window lenght for recursive exponential statistics
    # corresponding to 5% of the exponential maximum
    window = 0.5  # seconds
    samples = window / st[0].stats.delta
    # Compute the accumulation parameter lambda_
    lambda_ = 1 - math.exp(math.log(0.05) / samples)
    print(lambda_)
    lambda_ = 1. / samples
    print(lambda_)
    data_P, data_S, pol_filter =\
        rosenberger(st[2].data, st[1].data, st[0].data, lambda_)

    fig = plt.figure()
    ax1 = fig.add_subplot(411)
    ax1.set_title('Rosenberger C')
    ax2 = fig.add_subplot(412, sharex=ax1)
    ax3 = fig.add_subplot(413, sharex=ax1)
    ax4 = fig.add_subplot(414, sharex=ax1)
    ax4.set_xlabel('time (s)')

    ax1.set_ylim((-maxval, maxval))
    ax2.set_ylim((-maxval, maxval))
    ax3.set_ylim((-maxval, maxval))
    ax4.set_ylim((0, 1.2))

    ax1.plot(time, st[2].data, color='gray')
    ax1.plot(time, data_P[1], color='blue')
    ax1.plot(time, data_S[1], color='red')
    ax1.legend([st[2].stats.channel, 'P', 'S'])

    ax2.plot(time, st[1].data, color='gray')
    ax2.plot(time, data_P[2], color='blue')
    ax2.plot(time, data_S[2], color='red')
    ax2.legend([st[1].stats.channel, 'P', 'S'])

    ax3.plot(time, st[0].data, color='gray')
    ax3.plot(time, data_P[0], color='blue')
    ax3.plot(time, data_S[0], color='red')
    ax3.legend([st[0].stats.channel, 'P', 'S'])

    ax4.plot(time, pol_filter)

    plt.show()


if __name__ == '__main__':
    main()
