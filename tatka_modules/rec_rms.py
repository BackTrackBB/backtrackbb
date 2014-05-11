import os
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rec_rms.so')
lib_rec_rms = ctypes.CDLL(libpath)

lib_rec_rms._recursive_rms.argtypes = [
        ndpointer(ctypes.c_double),
        ndpointer(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_float,
        ]
lib_rec_rms._recursive_rms.restype = ctypes.c_void_p


def recursive_rms(signal, C_WIN):
    rms_signal = np.zeros(len(signal))
    lib_rec_rms._recursive_rms(
            signal, rms_signal, signal.size, C_WIN)
    return rms_signal


if __name__ == '__main__':
    signal = np.ones(60000)
    rms_signal = recursive_rms(signal, 1)
    print rms_signal
