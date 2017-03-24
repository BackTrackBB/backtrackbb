# -*- coding: utf8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from scipy.signal import lfilter
from math import floor, ceil
import numpy as np
import scipy as sp


def __gausscoeff(s):
    """Python implementation of the algorithm by Young & Vliet, 1995."""
    if s < .5:
        raise ValueError('Sigma for Gaussian filter must be >0.5 samples')
    q = 0.98711*s - 0.96330 if s > 0.5 else 3.97156 \
        - 4.14554*np.sqrt(1.0 - 0.26891*s)
    b = np.zeros(4)
    b[0] = 1.57825 + 2.44413*q + 1.4281*q**2 + 0.422205*q**3
    b[1] = 2.44413*q + 2.85619*q**2 + 1.26661*q**3
    b[2] = -(1.4281*q**2 + 1.26661*q**3)
    b[3] = 0.422205*q**3
    B = 1.0 - ((b[1] + b[2] + b[3])/b[0])
    # convert to a format compatible with lfilter's
    # difference equation
    B = np.array([B])
    A = np.ones(4)
    A[1:] = -b[1:]/b[0]
    return B, A


def Gaussian1D(signal, sigma, padding=0):
    n = signal.shape[0]
    tmp = np.zeros(n + padding)
    if tmp.shape[0] < 4:
        raise ValueError('Signal and padding too short')
    tmp[:n] = signal
    B, A = __gausscoeff(sigma)
    tmp = lfilter(B, A, tmp)
    tmp = tmp[::-1]
    tmp = lfilter(B, A, tmp)
    tmp = tmp[::-1]
    return tmp[:n]


def Gaussian2D(image, sigma, padding=0):
    n, m = image.shape[0], image.shape[1]
    tmp = np.zeros((n + padding, m + padding))
    if tmp.shape[0] < 4:
        raise ValueError('Image and padding too small')
    if tmp.shape[1] < 4:
        raise ValueError('Image and padding too small')
    B, A = __gausscoeff(sigma)
    tmp[:n, :m] = image
    tmp = lfilter(B, A, tmp, axis=0)
    tmp = np.flipud(tmp)
    tmp = lfilter(B, A, tmp, axis=0)
    tmp = np.flipud(tmp)
    tmp = lfilter(B, A, tmp, axis=1)
    tmp = np.fliplr(tmp)
    tmp = lfilter(B, A, tmp, axis=1)
    tmp = np.fliplr(tmp)
    return tmp[:n, :m]
#-----------------------------------------------------------------------------


def __shift(array, lag):
    pad = np.zeros(abs(lag))
    if lag > 0:
        return np.concatenate((pad, array[:-lag]))
    else:
        return np.concatenate((array[-lag:], pad))


def __shift2(array, lag):
    if lag == 0:
        return array
    ret = np.zeros_like(array)
    if lag > 0:
        ret[lag:] = array[:-lag]
    else:
        ret[:lag] = array[-lag:]
    return ret


def local_CC(sig1, sig2, t_lag, fs):
    """
    Calculation of Local-CC after Hale 2006.

    Output = H3 - non-smoothed LCC
          C - smoothed LCC (after Gussian smoothing)
          t_lag - array of t_lag times
    """
    l_max = int(t_lag*fs)
    h3 = np.zeros((2*l_max, len(sig1)), sig1.dtype)
    c = np.zeros((2*l_max, len(sig1)), sig1.dtype)
##
    t_lag = sp.linspace(-l_max/fs, l_max/fs, 2*l_max, endpoint=False)
##
##
    for l in xrange(-l_max, l_max):
        l_f = int(floor(l/2.))
        l_g = int(ceil(l/2.))
        h3[l+l_max] = (__shift2(sig1, l_f) * __shift2(sig2, -l_g) +
                       __shift2(sig1, l_g) * __shift2(sig2, -l_f))/2
        c[l+l_max] = Gaussian1D(h3[l+l_max], 10, padding=0)

#------------------not sure which formula for h3 is correct--------------------
#             h3[l+l_max]=(shift2(sig1,-l_f)*shift2(sig2,l_g) +\
#                    shift2(sig1,-l_g)*shift2(sig2,l_f))/2
#
    return c, h3, t_lag


def local_CCr(sig1, sig2, t_lag, fs, sigma):
    # Compute sigma in samples
    sigma = int(sigma * fs)
    l_max = int(t_lag*fs)
    h3 = np.zeros((2*l_max, len(sig1)), sig1.dtype)
    c = np.zeros((2*l_max, len(sig1)), sig1.dtype)
####
##
    for l in xrange(-l_max, l_max):
        l_f = int(floor(l/2.))
        l_g = int(ceil(l/2.))
        h3[l+l_max] = (__shift2(sig1, l_f)*__shift2(sig2, -l_g) +
                       __shift2(sig1, l_g)*__shift2(sig2, -l_f))/2
        c[l+l_max] = Gaussian1D(h3[l+l_max], sigma, padding=0)
#------------------check wich formula for h3 is correct--------------------
#             h3[l+l_max]=(shift2(sig1,-l_f)*shift2(sig2,l_g) +\
#                    shift2(sig1,-l_g)*shift2(sig2,l_f))/2
#
    return c, h3
