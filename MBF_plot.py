#!/usr/bin/env python
import sys
import os
import numpy as np
from tatka_modules.mod_filter_picker import make_LogFq,make_LinFq,MBfilter_CF,GaussConv
from tatka_modules.mod_setup import configure
from tatka_modules.read_traces import read_traces
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
from matplotlib.collections import LineCollection
import matplotlib.transforms as transforms

def main():
    line_width = 0.5

    config = configure()
    #---Reading data---------------------------------------------------------
    st_n = read_traces(config)
    st = st_n.select(station=config.stations[0])
    st.detrend(type='constant')
    st.detrend(type='linear')
    print st
    hos_sigma = config['hos_sigma_' + 'S']
    if config.cut_data:
        st.trim(st[0].stats.starttime + config.cut_start,
                st[0].stats.starttime + config.cut_start + config.cut_delta)

    dt1 = st[0].stats.delta

    if config.band_spacing == 'lin':
        frequencies = make_LinFq(config.f_min, config.f_max, dt1, config.n_freq_bands)
    elif config.band_spacing == 'log':
        frequencies = make_LogFq(config.f_min, config.f_max, dt1, config.n_freq_bands)
    frequencies = frequencies[::-1]

    #---------------------------MB filtering------------------------------------
    HP2_1, MBkurt_1, Tn2_1, Nb2_1 = MBfilter_CF(st, frequencies,
                                         var_w=config.win_type,
                                         CF_type=config.ch_function,
                                         order=config.hos_order,
                                         CF_decay_win=config.decay_const,
                                         filter_strength=config.filter_strength,
                                         hos_sigma=hos_sigma[config.stations[0]])
    print 'Creating characteristic function: %s' % (st[0].stats.station)
    if config.ch_function == 'kurtosis':
        MBkurt_1max_gauss = GaussConv(np.amax(MBkurt_1,axis=0),int(config.decay_const/dt1/2))
        MBkurt_1max = np.amax(MBkurt_1,axis=0)
    elif config.ch_function == 'envelope' or config.ch_function == 'hilbert':
        MBkurt_1max = np.sqrt(np.power(MBkurt_1, 2).mean(axis=0))
    ###---------------------------------------------------------------------------

    #-------                  Plotting
    fig = plt.figure(figsize=(18,10))
    # ax1: original trace
    ax1 = fig.add_axes([0.02, 0.207, 0.47, 0.17])
    # ax2: filtered traces
    ax2 = fig.add_axes([0.02, 0.395, 0.47, 0.6], sharex=ax1)
    # ax3: characteristic functions
    ax3 = fig.add_axes([0.51, 0.395, 0.47, 0.6], sharex=ax1, sharey=ax2)
    # ax4: summary characteristic function
    ax4 = fig.add_axes([0.51, 0.207, 0.47, 0.17], sharex=ax1, sharey=ax1)

    time_array = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate

    ax1.set_ylim(-1, 1)
    ax4.set_ylim(-1, 1)
    ax2.set_xlim(0, max(time_array))
    ax3.set_xlim(0, max(time_array))
    ax2.set_ylim(-1, Nb2_1)
    ax3.set_ylim(-1, Nb2_1)
    ax1.set_xlabel('Time [s]',fontsize=14)
    ax4.set_xlabel('Time [s]',fontsize=14)
    ax1.yaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax3.yaxis.set_visible(False)
    ax2.xaxis.set_visible(False)
    ax3.xaxis.set_visible(False)
    ax4.yaxis.set_visible(False)


    # ax1: original trace:
    st.normalize()
    if len(st) == 2:
        ax1.plot(time_array, st[0].data,'k',lw=line_width,label = st[0].id)
        ax1.plot(time_array, st[1].data-max(st[0].data),'b',lw=line_width,
                 label = st[1].id)
    else:
        ax1.plot(time_array, st[0].data,'k',lw=line_width,label = st[0].id)
    ax1.legend()
    # end ax1


    # ax2, ax3: filtered traces, characteristic functions
    # Normalize to one
    HP2_1 /= np.amax(HP2_1)
    MBkurt_1 /= np.amax(MBkurt_1)

    # Create line objects and add to plot
    traces = [trace + i for i, trace in enumerate(list(HP2_1))]
    lines = LineCollection([list(zip(time_array, trace)) for trace in traces],
            linewidth=line_width, color='k', rasterized=True)
    ax2.add_collection(lines)

    traces = [trace + i for i, trace in enumerate(list(MBkurt_1))]
    lines = LineCollection([list(zip(time_array, trace)) for trace in traces],
            linewidth=line_width, color='k', rasterized=True)
    ax3.add_collection(lines)

    # Transformations for label positioning
    trans1 = transforms.blended_transform_factory(
        fig.transFigure, ax2.transData)
    trans2 = transforms.blended_transform_factory(
        fig.transFigure, ax3.transData)
    inv = fig.transFigure.inverted()
    xtext1 = inv.transform(ax2.transData.transform((time_array[0], 0)))[0]
    xtext2 = inv.transform(ax3.transData.transform((time_array[0], 0)))[0]

    for i, Tn in enumerate(Tn2_1):
        label = '%.2f Hz' % (1./Tn)
        ax2.text(xtext1, i+0.1, label, fontsize=13, transform=trans1, clip_on=True)
        ax3.text(xtext2, i+0.1, label, fontsize=13, transform=trans2, clip_on=True)
    # end ax2, ax3


    # ax4: summary characteristic function
    label = 'MBF sum. ' + config.ch_function
    ax4.plot(time_array, MBkurt_1max/max(MBkurt_1max),'r', label=label, lw=1.5)
    if config.ch_function == 'kurtosis':
        label = 'MBF sum. gauss'
        ax4.plot(time_array, MBkurt_1max_gauss/max(MBkurt_1max_gauss), 'b', label=label, lw=1.5)

    st.filter('bandpass', freqmin=config.f_min,
              freqmax=config.f_max, corners=2, zerophase=False)
    st.normalize()
    if len(st) == 2:
        horiz_env = np.sqrt(np.power(st[0].data, 2) +
                            np.power(st[1].data, 2))
        ax4.plot(time_array, horiz_env,'k',alpha = 0.6,lw=0.5)
    else:
        ax4.plot(time_array, st[0].data,'k',alpha = 0.6,lw=0.5)
    ax4.legend()
    # end ax4
    #------------------------------------------------------------------------------

    if len(st) == 2:
        name_fig_out = 'MBFplot.' + config.ch_function + '.' + st[0].stats.station + '.horiz'
    if len(st) == 1:
        name_fig_out = 'MBFplot.' + config.ch_function + '.' + st[0].id

    #fig.savefig(os.path.join(config.out_dir,name_fig_out), format=config.plot_format)
    plt.show()


if __name__ == '__main__':
    main()
