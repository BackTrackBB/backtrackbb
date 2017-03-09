# -*- coding: utf8 -*-
# Wavefield separation (P and non-P) based on polarization analysis.
#
# Rosenberger, A., 2010. Real-Time Ground-Motion Analysis:
# Distinguishing P and S Arrivals in a Noisy Environment.
# Bull. Seismol. Soc. Am. 100, 1252–1262. doi: 10.1785/0120090265
#
# (c) 2014 - Pierre Romanet <romanet@ipgp.fr>,
#            Claudio Satriano <satriano@ipgp.fr>
import numpy as np


def _update_(U, D, d, lambda_):
    """Go from u(n) to u(n+1)."""
    I = np.identity(3)

    m = U.T.dot(d)
    p = (I - U.dot(U.T)).dot(d)
    p_norm = np.linalg.norm(p)

    # Make p and m column vectors
    p = p[np.newaxis].T
    m = m[np.newaxis].T

    U_left = np.hstack((U, p/p_norm))
    Q = np.hstack((lambda_ * D, m))
    Q = np.vstack((Q, [0, 0, p_norm]))

    # SVD
    U_right, D_new, V_left = np.linalg.svd(Q)

    # Get rid of the smallest eigenvalue
    D_new = D_new[0:2]
    D_new = np.diagflat(D_new)

    U_right = U_right[:, 0:2]

    return U_left.dot(U_right), D_new


def rosenberger(dataX, dataY, dataZ, lambda_):
    """
    Separate P and non-P wavefield from 3-component data.

    Return a two set of 3-component traces.
    """
    # Construct the data matrix
    A = np.vstack((dataZ, dataX, dataY))

    # SVD of the first 3 samples:
    U, D, V = np.linalg.svd(A[:, 0:3])

    # Get rid of the smallest eigenvalue
    D = D[0:2]
    D = np.diagflat(D)
    U = U[:, 0:2]

    save_U = np.zeros(len(dataX))
    save_U[0] = abs(U[0, 0])

    Dp = np.zeros((3, len(dataX)))
    Ds = np.zeros((3, len(dataX)))
    Dp[:, 0] = abs(U[0, 0]) * A[:, 2]
    Ds[:, 0] = (1 - abs(U[0, 0])) * A[:, 2]

    # Loop over all the values
    for i in range(1, A.shape[1]):
        d = A[:, i]
        U, D = _update_(U, D, d, lambda_)

        Dp[:, i] = abs(U[0, 0]) * d
        Ds[:, i] = (1-abs(U[0, 0])) * d

        save_U[i] = abs(U[0, 0])

    return Dp, Ds, save_U


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
    window = 0.5  # seconds
    samples = window / st[0].stats.delta
    # Compute the accumulation parameter lambda_
    lambda_ = math.exp(math.log(0.05) / samples)
    data_P, data_S, U =\
        rosenberger(st[2].data, st[1].data, st[0].data, lambda_)

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
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
