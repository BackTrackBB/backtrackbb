import numpy as np
import scipy as sp
from obspy.signal.util import smooth
#from obspy.signal.invsim import cosTaper
from rec_filter import recursive_filter
from rec_rms import recursive_rms
from rec_kurtosis import recursive_kurtosis
from RosenbergerAlgorithm import rosenberger
from scipy.signal import gaussian
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


def MBfilter_CF(st, fq, n_win, CF_type='envelope', var_w=True, 
                C_kurtosis=0.01, order1=4, order2=2, power2=2,
                C_rosenberger=0.99):
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

# Less than 3 components------------------------------------------------------------------------
    if len(st) < 3:
        
        y[0]=np.mean(y[0:n_win])

        #taper=cosTaper(len(y),p=0.1)
        #y=y*taper

        # derivative
        dy = np.gradient(y)

        YN1 = np.zeros((Nb, len(y)), float)
        CF = np.zeros((Nb, len(y)), float)

        for n in xrange(Nb):
            YN1[n] = recursive_filter(dy, CN_HP[n], CN_LP[n])

            if var_w:
                n_win_mb=(b*Tn[n]/dT)
                #print n_win_mb, n_win
            else:
                n_win_mb = (n_win)
    
            if CF_type == 'envelope':
                #module using python functions
                ##CF[n] = meansq_rec(YN1[n],n_win_mb)

                #module using C function
                CF[n] = recursive_rms(YN1[n],1./n_win_mb)

            if CF_type == 'hilbert':
                CF[n] = smooth(abs(sp.signal.hilbert(YN1[n])),n_win_mb)

            if CF_type == 'kurtosis':
                #module using python functions
                ##CF[n] = recKurt_1(YN1[n], n_win_mb)

                #module using C function
                CF[n] = recursive_kurtosis(YN1[n], C_kurtosis, 0.001, 4, 2, 2) 

# More than 3 components------------------------------------------------------------------------
    else:
        
        tr2 = st.select(component='E')[0]
        tr3 = st.select(component='N')[0]

        y2 = tr2.data
        y3 = tr3.data

        y2 = y2 - y2.mean()
        y2[0] = np.mean(y2[0:n_win])

        y3 = y3 - y3.mean()
        y3[0] = np.mean(y3[0:n_win])

        #taper=cosTaper(len(y),p=0.1)
        #y=y*taper

        # derivative
        dy1 = np.gradient(y)
        dy2 = np.gradient(y2)
        dy3 = np.gradient(y3)

        # Initializing arrays
        YN1 = np.zeros((Nb, len(y)), float)
        YN2 = np.zeros((Nb, len(y)), float)
        YN3 = np.zeros((Nb, len(y)), float)
        CF = np.zeros((Nb, len(y)), float)

        for n in xrange(Nb):
            
            YN1[n] = recursive_filter(dy1, CN_HP[n], CN_LP[n])
            YN2[n] = recursive_filter(dy2, CN_HP[n], CN_LP[n])
            YN3[n] = recursive_filter(dy3, CN_HP[n], CN_LP[n])

            print 'Rosenberger in process {}/{}'.format(n+1,Nb)
            filtered_dataP, filtered_dataS, U = rosenberger(YN2[n], YN3[n], YN1[n], C_rosenberger)
            

            if var_w:
                n_win_mb = (b*Tn[n]/dT)
                #print n_win_mb, n_win
            else:
                n_win_mb = (n_win)

            if CF_type == 'kurtosis':
                #module using python functions
                ##CF[n] = recKurt_1(YN[n], n_win_mb)

                #module using C function
                CF[n,2:] = recursive_kurtosis(filtered_dataP[0,:], C_kurtosis, 0.001, 4, 2, 2) 

#----------------------------------------------------------
    return YN1, CF, Tn, Nb

def GaussConv(data_in,sigma):
    derivative_data = np.zeros(len(data_in),float)
    for i in xrange(0,len(data_in)-1):
        derivative_data[i] = data_in[i+1]-data_in[i]
    
    derivative_data [derivative_data< 0] = 0
    gauss_window = gaussian(len(derivative_data),sigma)
    CF_gaussian = np.convolve(derivative_data, gauss_window,mode='same')
    
    return CF_gaussian



if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    from generate_signal import *
    from obspy.core import read, Trace, Stream

    #if arguments, read the file
    if len(sys.argv) == 2:
        filename = sys.argv[1]
        print filename
        data = read(filename)
        signal = np.array(data[0].data,dtype='float')
        sys.exitfunc()
    elif len(sys.argv) >= 2:
        #to be completed to account for extra arguments
        pass
    else:
        noise = generate_signal_noise2(1000, 0.05)
        signal = generate_signal_expSin(300,0.005,0.5,noise,0.5,500,0.05,1)

    # sampling frequency
    if 'data' in locals():
        delta = data[0].stats.sampling_rate
    else:
        delta = 0.01
    print 'Sampling frequency = {}'.format(delta)

    tr = Trace(signal)
    tr.stats.delta = delta
    tr.stats.channel = 'HHZ'
    st = Stream(tr)
    
    fq = make_LinFq(1, 40, delta, 8)
    YN, CF, Tn, Nb = MBfilter_CF(st, fq, 20, CF_type='kurtosis', var_w=False)

    fig1 = plt.figure()
    fig2 = plt.figure()
    max2 = CF.max()

    for n in range(Nb):
        ax1 = fig1.add_subplot(Nb+1, 1, n+1)
        ax2 = fig2.add_subplot(Nb+1, 1, n+1)
	    #ax2.set_ylim((0, max2))
        ax1.plot(YN[n])
        ax2.plot(CF[n])

    ax1.plot(signal,'g')
    ax2.plot(recursive_kurtosis(signal, 0.1, 0.001, 4, 2, 2),'g')

    plt.show()
