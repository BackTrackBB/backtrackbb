import os
from ctypes import CDLL, c_int, c_float, c_double, c_void_p, POINTER, byref
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rec_filter.so')
lib_rec_filter = CDLL(libpath)

lib_rec_filter._recursive_filter_BP.argtypes = [
        ndpointer(dtype=np.float64), #signal
        ndpointer(dtype=np.float64), #filt_signal
        c_int, #npts
        c_float, #C_HP
        c_float, #C_LP
        POINTER(c_double), #filterH1
        POINTER(c_double), #filterH2
        POINTER(c_double), #filterL1
        POINTER(c_double) #filterL2
        ]
lib_rec_filter._recursive_filter_BP.restype = c_void_p

lib_rec_filter._recursive_filter_HP.argtypes = [
        ndpointer(dtype=np.float64), #signal
        ndpointer(dtype=np.float64), #filt_signal
        c_int, #npts
        c_float, #C_HP
        POINTER(c_double), #filterH1
        POINTER(c_double) #filterH2
        ]
lib_rec_filter._recursive_filter_HP.restype = c_void_p


def recursive_filter(signal, C_HP, C_LP=None, rec_memory=None):
    signal = np.array(signal, dtype=np.float64)
    filt_signal = np.zeros(len(signal))

    if rec_memory is not None:
        filterH1 = c_double(rec_memory.filterH1)
        filterH2 = c_double(rec_memory.filterH2)
        filterL1 = c_double(rec_memory.filterL1)
        filterL2 = c_double(rec_memory.filterL2)
    else:
        filterH1 = c_double(0)
        filterH2 = c_double(0)
        filterL1 = c_double(0)
        filterL2 = c_double(0)

    if C_LP is not None:
        lib_rec_filter._recursive_filter_BP(
            signal, filt_signal, signal.size, C_HP, C_LP,
            byref(filterH1), byref(filterH2),
            byref(filterL1), byref(filterL2))
    else:
        lib_rec_filter._recursive_filter_HP(
            signal, filt_signal, signal.size, C_HP,
            byref(filterH1), byref(filterH2))

    if rec_memory is not None:
        rec_memory.filterH1 = filterH1.value
        rec_memory.filterH2 = filterH2.value
        rec_memory.filterL1 = filterL1.value
        rec_memory.filterL2 = filterL2.value

    return filt_signal


if __name__ == '__main__':
    signal = np.ones(60000)
    filt_signal = recursive_filter(signal, 0.1, 0.1)
    print filt_signal
