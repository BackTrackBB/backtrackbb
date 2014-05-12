import os
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np


libpath = os.path.join(os.path.dirname(__file__), os.pardir, 'lib', 'lib_rec_cc.so')
lib_rec_cc = ctypes.CDLL(libpath)

lib_rec_cc._local_CCr.argtypes = [
        ndpointer(dtype=np.float64),
        ndpointer(dtype=np.float64),
        ctypes.c_int,
        ndpointer(dtype=np.float64),
        ndpointer(dtype=np.float64),
        ctypes.c_int,
        ctypes.c_double
        ]
lib_rec_cc._local_CCr.restype = ctypes.c_void_p


def local_CCr(signal1, signal2, t_lag, fs, sigma):
    signal1 = np.array(signal1, dtype=np.float64)
    signal2 = np.array(signal2, dtype=np.float64)
    if signal1.size != signal2.size:
        raise RuntimeError, 'Signals must have the same size.'

    # Compute sigma in samples
    sigma *= fs
    if sigma < 0.5:
        raise RuntimeError, 'Sigma for Gaussian filter must be >=0.5 samples.'
    lmax = int(t_lag * fs)
    cc_no_filt = np.zeros(2*lmax * len(signal1))
    cc = np.zeros(2*lmax * len(signal1))
    lib_rec_cc._local_CCr(
            signal1, signal2, signal1.size,
            cc_no_filt, cc, lmax,
            sigma)
    return cc.reshape(2*lmax, len(signal1)), cc_no_filt.reshape(2*lmax, len(signal1))


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    signal1 = np.zeros(6001)
    signal2 = np.zeros(6001)
    signal1[1000:1200] = 1
    signal2[1200:1300] = 1
    cc, cc_no_filt = local_CCr(signal1, signal2, 300, 1, 100)
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.imshow(cc_no_filt)
    ax2 = fig.add_subplot(212, sharex=ax1, sharey=ax1)
    ax2.imshow(cc)
    plt.show()
