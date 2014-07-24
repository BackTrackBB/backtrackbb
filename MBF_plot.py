#!/usr/bin/env python
import sys
import os
import numpy as np
from obspy.core import UTCDateTime, read
from tatka_modules.mod_filter_picker import make_LogFq,make_LinFq,MBfilter_CF,GaussConv
from tatka_modules.parse_config import parse_config
from tatka_modules.read_traces import read_traces
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
matplotlib.rcParams['pdf.fonttype'] = 42

line_width = 0.5

if len(sys.argv) != 2:
    print "this_code  <input_config_file>"
    sys.exit(1)
else:
    Config_file = sys.argv[1]
if not os.path.isfile(Config_file):
    print "File {0} does not exist".format(Config_file)
    sys.exit(1)

config = parse_config(Config_file)
st_n, station = read_traces(config)

st = st_n.select(station=config.stations[0])     
st.detrend(type='constant')
st.detrend(type='linear')
print st
hos_sigma = config['hos_sigma_' + 'S']
if config.cut_data:
    st.trim(st[0].stats.starttime+config.cut_start,
            st[0].stats.starttime+config.cut_start+config.cut_delta)

ttime = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate
dt1=st[0].stats.delta

if config.band_spacing == 'lin':
    frequencies = make_LinFq(config.f_min, config.f_max, dt1, config.n_freq_bands)
elif config.band_spacing == 'log':
    frequencies = make_LogFq(config.f_min, config.f_max, dt1, config.n_freq_bands)
frequencies = frequencies[::-1]

#---------------------------MB filtering------------------------------------
HP2_1,  MBkurt_1, Tn2_1, Nb2_1 = MBfilter_CF(st, frequencies,
                                     var_w = config.varWin_stationPair,
                                     CF_type = config.ch_function,
                                     CF_decay_win = config.decay_const,
                                     hos_sigma=hos_sigma[config.stations[0]])
print 'Creating characteristic function: %s' % (st[0].stats.station)
if config.ch_function == 'kurtosis':         
    MBkurt_1max = GaussConv(np.amax(MBkurt_1,axis=0),int(config.decay_const/dt1/4))
elif config.ch_function == 'envelope' or config.ch_function == 'hilbert': 
    MBkurt_1max = np.sqrt(np.power(MBkurt_1, 2).mean(axis=0))
###---------------------------------------------------------------------------
#-------                  Plotting
fig = plt.figure(figsize=(18,10))
ax1 = fig.add_axes([0.02, 0.395, 0.47, 0.6])
ax2 = fig.add_axes([0.51, 0.395, 0.47, 0.6],sharex=ax1)
ax3 = fig.add_axes([0.515, 0.207, 0.47, 0.17],sharex=ax1)
ax4 = fig.add_axes([0.02, 0.207, 0.47, 0.17],sharex=ax1)

ax1.set_xlim(0,max(ttime))
ax2.set_xlim(0,max(ttime))
ax1.yaxis.set_visible(False)
ax2.yaxis.set_visible(False)
ax1.xaxis.set_visible(False)
ax2.xaxis.set_visible(False)
ax3.yaxis.set_visible(False)
ax4.yaxis.set_visible(False)
    
scale1=np.amax(HP2_1)
scale2=np.amax(MBkurt_1)

k=0
for i in xrange(0,len(Tn2_1)):
    k +=1
    ax1.plot(ttime,HP2_1[i,:]+scale1*(k-1),'k',label="CF",lw=line_width)  
    ax2.plot(ttime,MBkurt_1[i,:]+scale2*(k-1),'k',lw=line_width)
    note=str('%.2f' % (1./Tn2_1[i])+' Hz')
    ax1.text(min(ttime),scale1*(k-1)-scale1/3, note,fontsize=13)
    ax2.text(min(ttime),scale2*(k-1)-scale2/3, note,fontsize=13)


sta1 = st[0].stats.station+'.'+st[0].stats.channel
s_rate1 =str(st[0].stats.sampling_rate)

y_min1=-abs(max(HP2_1[0]))
y_min2=-abs(max(MBkurt_1[0]))
y_max1=scale1*(Nb2_1)
y_max2=scale2*(Nb2_1)
ax1.set_ylim(y_min1, y_max1)
ax2.set_ylim(y_min2, y_max2)

ax3.plot(ttime, MBkurt_1max/max(MBkurt_1max),'r',label="MBF Sum. "+config.ch_function,lw=1.5)
if len(st)==2:
    ax4.plot(ttime, st[0].data,'k',lw=line_width,label = st[0].id)
    ax4.plot(ttime, st[1].data-max(st[0].data),'b',lw=line_width,
             label = st[1].id)
else:
    ax4.plot(ttime, st[0].data,'k',lw=line_width,label = st[0].id)

st.filter('bandpass', freqmin=config.f_min,
          freqmax=config.f_max, corners=2, zerophase=False)
st.normalize()
if len(st)==2:
    horiz_env = np.sqrt(np.power(st[0].data,2)+
                        np.power(st[1].data,2))
    ax3.plot(ttime, horiz_env,'k',alpha = 0.6,lw=0.5)
else:
    ax3.plot(ttime, st[0].data,'k',alpha = 0.6,lw=0.5)    
#------------------------------------------------------------------------------
ax3.set_ylim(-1, 1)
ax3.set_xlabel('Time [s]',fontsize=14)
ax4.set_xlabel('Time [s]',fontsize=14)
ax3.legend()
ax4.legend()

if len(st) == 2:
    name_fig_out = 'MBFplot.'+config.ch_function+'.'+st[0].stats.station+'.horiz'
if len(st) == 1:
    name_fig_out = 'MBFplot.'+config.ch_function+'.'+st[0].id
    
fig.savefig(os.path.join(config.out_dir,name_fig_out), format=config.plot_format)
#plt.show()

