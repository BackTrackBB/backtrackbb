import os
from ctypes import CDLL, c_int, c_float, c_double, c_void_p, POINTER, byref
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rec_rms.so')
lib_rec_rms = CDLL(libpath)

lib_rec_rms._recursive_rms.argtypes = [
        ndpointer(dtype=np.float64), #signal
        ndpointer(dtype=np.float64), #rms_signal
        c_int, #npts
        c_float, #C_WIN
        POINTER(c_double) #mean_sq
        ]
lib_rec_rms._recursive_rms.restype = c_void_p


def recursive_rms(signal, C_WIN, rec_memory=None):
    signal = np.array(signal, dtype=np.float64)
    rms_signal = np.zeros(len(signal))

    if rec_memory is not None:
        mean_sq = c_double(rec_memory.mean_sq)
    else:
        mean_sq = c_double(0)

    lib_rec_rms._recursive_rms(
            signal, rms_signal, signal.size, C_WIN, byref(mean_sq))

    if rec_memory is not None:
        rec_memory.mean_sq = mean_sq.value

    return rms_signal


if __name__ == '__main__':
    signal = np.ones(60000)
    rms_signal = recursive_rms(signal, 1)
    print rms_signal
