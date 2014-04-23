import numpy as np
import scipy as sp
from obspy.signal.util import smooth
#from obspy.signal.invsim import cosTaper
from rec_filter import recursive_filter
from rec_rms import recursive_rms
from rec_kurtosis import recursive_kurtosis
#------------------------------------------------------------------------------------------------------------


def make_LinFq(f_min, f_max, delta, nfreq):
    """
    Calculates linearly spaced frequency array for MB filtering
    """
    f_ny = float(1./(2*delta))
    if f_max > f_ny:
        f_max = f_ny

    freq = np.linspace(f_min, f_max, nfreq)
    return freq[::-1]



def make_LogFq(f_min, f_max, delta, nfreq):
    """
    Calculates log spaced frequency array for MB filtering
    """
    f_ny = float(1./(2*delta))
    if f_max > f_ny:
        f_max = f_ny

    freq = np.logspace(np.log2(f_min), np.log2(f_max), nfreq, base=2)
    return freq[::-1]


def MBfilter_CF(st, fq, n_win, CF_type='envelope', var_w=True):
    """
    Performs MBfiltering using 2HP+2LP recursive filter
    and calculates the characteristic function (CF)
    for each band.
    """
    tr = st.select(component='Z')[0]
    y = tr.data
    dT = tr.stats.delta
    Nb = len(fq)
    Tn = 1/fq
    wn = Tn/(2*np.pi)         #wn=Tn/(2*pi)- time constant
    CN_HP = wn/(wn+dT)        # hight-pass filter constant
    CN_LP = dT/(wn+dT)        # low-pass filter constant

    b = n_win*dT/Tn[int(Nb/2+1)]

    y = y - y.mean()
    y[0]=np.mean(y[0:n_win])

    #taper=cosTaper(len(y),p=0.1)
    #y=y*taper

    # derivative
    dy = np.gradient(y)

    YN=np.zeros((Nb, len(y)), float)
    CF=np.zeros((Nb, len(y)), float)

    for n in xrange(Nb):
        YN[n] = recursive_filter(dy, CN_HP[n], CN_LP[n])

        if var_w:
            n_win_mb=(b*Tn[n]/dT)
                #print n_win_mb, n_win
        else:
            n_win_mb = (n_win)

        if CF_type == 'envelope':
            #module using python functions
            ##CF[n] = meansq_rec(YN[n],n_win_mb)

            #module using C function
            CF[n] = recursive_rms(YN[n],1./n_win_mb)

        if CF_type == 'hilbert':
            CF[n] = smooth(abs(sp.signal.hilbert(YN[n])),n_win_mb)

        if CF_type == 'kurtosis':
            #module using python functions
            ##CF[n] = recKurt_1(YN[n], n_win_mb)

            #module using C function
            CF[n] = recursive_kurtosis(YN[n], 1./n_win_mb, 0.01, 4, 2, 2)

#----------------------------------------------------------
    return YN, CF, Tn, Nb
