# -*- coding: utf8 -*-
# rosenberger.py
#
# Python interface to C code for P and non-P wavefield separation
# by Andreas Rosenberger.
#
# Rosenberger, A., 2010. Real-Time Ground-Motion Analysis: Distinguishing P and S Arrivals
# in a Noisy Environment. Bull. Seismol. Soc. Am. 100, 1252–1262. doi: 10.1785/0120090265
#
# (c) 2014 - Claudio Satriano <satriano@ipgp.fr>, Pierre Romanet <romanet@ipgp.fr>
import os
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rosenberger.so')
lib_rosenberger = ctypes.CDLL(libpath)

lib_rosenberger.rosenberger.argtypes = [
        ndpointer(dtype=np.float64), #dataX
        ndpointer(dtype=np.float64), #dataY
        ndpointer(dtype=np.float64), #dataZ
        ndpointer(dtype=np.float64), #dataX_P
        ndpointer(dtype=np.float64), #dataY_P
        ndpointer(dtype=np.float64), #dataZ_P
        ndpointer(dtype=np.float64), #dataX_S
        ndpointer(dtype=np.float64), #dataY_S
        ndpointer(dtype=np.float64), #dataZ_S
        ctypes.c_int, #npts
        ctypes.c_float, #lambda
        ctypes.c_float, #delta
        ctypes.c_char, #proj
        ctypes.c_char #rl_filter
        ]
lib_rosenberger.rosenberger.restype = ctypes.c_void_p


def rosenberger(dataX, dataY, dataZ,
                lambda_, delta=0.0, proj=False, rl_filter=False):
    """
       Separates P and non-P wavefield from 3-component data
       and returns it as two set of 3-component traces.
    """
    dataX = np.array(dataX, dtype=np.float64)
    dataY = np.array(dataY, dtype=np.float64)
    dataZ = np.array(dataZ, dtype=np.float64)

    dataX_P = np.zeros_like(dataX)
    dataY_P = np.zeros_like(dataY)
    dataZ_P = np.zeros_like(dataZ)
    dataX_S = np.zeros_like(dataX)
    dataY_S = np.zeros_like(dataY)
    dataZ_S = np.zeros_like(dataZ)

    delta = ctypes.c_float(delta)
    proj = chr(int(proj))
    rl_filter = chr(int(rl_filter))

    lib_rosenberger.rosenberger(dataX, dataY, dataZ,
                                dataX_P, dataY_P, dataZ_P,
                                dataX_S, dataY_S, dataZ_S,
                                len(dataX),
                                lambda_, delta, proj, rl_filter)

    data_P = np.array([dataZ_P, dataX_P, dataY_P])
    data_S = np.array([dataZ_S, dataX_S, dataY_S])
    return data_P, data_S, None


def main():
    import math
    from obspy import read
    import matplotlib.pyplot as plt

    # We use the default ObsPy example
    st = read()
    st.filter(type='highpass', freq=1.0)
    maxval = max(st.max())

    time = np.arange(len(st[1].data)) * st[1].stats.delta

    # Window lenght for recursive exponential statistics
    # corresponding to 5% of the exponential maximum
    window = 0.5 #seconds
    samples = window / st[0].stats.delta
    # Compute the accumulation parameter lambda_
    lambda_ = math.exp(math.log(0.05) / samples)
    data_P, data_S, U = rosenberger(st[2].data, st[1].data, st[0].data, lambda_)

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.set_title('Rosenberger C')
    ax2 = fig.add_subplot(312, sharex=ax1)
    ax3 = fig.add_subplot(313, sharex=ax1)
    ax3.set_xlabel('time (s)')

    ax1.set_ylim((-maxval, maxval))
    ax2.set_ylim((-maxval, maxval))
    ax3.set_ylim((-maxval, maxval))

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

    plt.show()


if __name__ == '__main__':
    main()
