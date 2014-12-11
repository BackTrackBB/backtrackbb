import os
import numpy as np
import pylab
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.axes_grid1 import make_axes_locatable


def bp_plot(config, proj_grid,
            coord_eq, t_b, t_e, datestr, fq_str,
            coord_sta,
            st, sta, st_CF,
            time, time_env,
            fq, n1, n22, trigger,
            arrival_times=None, Mtau=None):

    LTrig = config.trigger
    lcc_max = config.lcc_max
    out_dir = config.out_dir
    scmap = config.scmap
    plot_waveforms = config.plot_waveforms

    Xmin, Xmax, Ymin, Ymax, Zmin, Zmax = proj_grid.get_extent()
    ratio = (Xmax - Xmin) / (Ymax - Ymin)
    fig = figure.Figure(figsize=(20, 20))
    plot_xz_size = ((Zmax - Zmin)/(Xmax - Xmin))*100
    plot_yz_size = plot_xz_size / ratio
    plot_cbar_size = 5 #percent
    xz_size = '%f %%' % plot_xz_size
    yz_size = '%f %%' % plot_yz_size
    cb_size = '%f %%' % plot_cbar_size
    sta_smbl_size = 150 / ratio
    eq_smbl_size = 200 / ratio
    trig_smbl_size = 200 / ratio

    scmap2 = pylab.cm.jet
    scmap2.set_under('w', LTrig)
    lcc_min = LTrig

    if trigger is not None:
        i_max, j_max, k_max = map(int, (trigger.i, trigger.j, trigger.k))
        xx_max, yy_max, zz_max = trigger.x, trigger.y, trigger.z

#--ax1: stacked grid
#--ax1_xy
    ax1_xy = fig.add_subplot(221)
    divider1 = make_axes_locatable(ax1_xy)

    if trigger is not None:
        i_grid = i_max
        j_grid = j_max
        k_grid = k_max
    elif coord_eq:
        x_eq, y_eq, z_eq = coord_eq
        i_grid, j_grid, k_grid = proj_grid.get_ijk(x_eq[0], y_eq[0], z_eq[0])
    else:
        i_grid = int(proj_grid.nx / 2.)
        j_grid = int(proj_grid.ny / 2.)
        k_grid = int(proj_grid.nz / 2.)

    xx_grid, yy_grid, zz_grid = proj_grid.get_xyz(i_grid, j_grid, k_grid)

    hnd=ax1_xy.imshow(np.flipud(np.transpose(proj_grid[:,:,k_grid])),
                             extent=proj_grid.get_xy_extent(), cmap=scmap, rasterized=True)
    ax1_xy.axis('tight')
    ax1_xy.set_xlim(Xmin,Xmax)
    ax1_xy.set_ylim(Ymin,Ymax)
    ax1_xy.set_xlabel('X[km]')
    ax1_xy.set_ylabel('Y[km]')
    labels = ax1_xy.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    tt1 = st[0].stats.starttime+t_b
    tt2 = st[0].stats.starttime+t_e
    ax1_xy.set_title('Date: %s, Time: %s.%03d - %s.%03d (%s - %s s), depth: %s km' %
                     (tt1.date,
                      tt1.strftime('%H:%M:%S'),
                      int(round(tt1.microsecond/1000.)),
                      tt2.strftime('%H:%M:%S'),
                      int(round(tt2.microsecond/1000.)),
                      t_b + config.cut_start,
                      t_e + config.cut_start,
                      zz_grid))
    if coord_eq:
        ax1_xy.scatter(coord_eq[0], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1,c='w')
    for sta in coord_sta:
        x_sta, y_sta = coord_sta[sta]
        ax1_xy.scatter(x_sta, y_sta, marker='^', s=sta_smbl_size, linewidths=1, c='k', alpha=0.79)
        trans = ax1_xy.transData + ax1_xy.transAxes.inverted()
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        ax1_xy.text(x_sta_ax+0.02, y_sta_ax+0.02, sta, fontsize=12, color='k', transform=ax1_xy.transAxes)

    if trigger is not None:
        ax1_xy.scatter(xx_max, yy_max,
                    marker='*', s=trig_smbl_size, linewidths=1,c='g')
    ax1_xy.set_aspect('equal', 'datalim')

#--ax1_yz
    ax1_yz = divider1.append_axes('right', size=yz_size, pad=0.05, sharey=ax1_xy)
    ax1_yz.imshow(np.flipud(proj_grid[i_grid,:,:]),
                     extent=proj_grid.get_zy_extent(), cmap=scmap, rasterized=True)

    if trigger is not None:
        ax1_yz.scatter(zz_max, yy_max,
                    marker='*', s=trig_smbl_size, linewidths=1, c='g')
    ax1_yz.axis('tight')
    ax1_yz.set_xlim(Zmin,Zmax)
    ax1_yz.set_ylim(Ymin,Ymax)
    ax1_yz.set_xlabel('Z[km]')
    labels = ax1_yz.get_xticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    labels = ax1_yz.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    ax1_yz.yaxis.set_visible(False)

#--ax1_xz
    ax1_xz = divider1.append_axes('bottom', size=xz_size, pad=0.05, sharex=ax1_xy)
    ax1_xz.imshow(np.flipud(np.transpose(proj_grid[:,j_grid,:])),
                     extent=proj_grid.get_xz_extent(), cmap=scmap, rasterized=True)

    if trigger is not None:
        ax1_xz.scatter(xx_max, zz_max,
                    marker='*', s=trig_smbl_size, linewidths=1, c='g')

    ax1_xz.axis('tight')
    ax1_xz.set_xlim(Xmin,Xmax)
    ax1_xz.set_ylim(Zmin,Zmax)
    ax1_xz.set_ylim(ax1_xz.get_ylim()[::-1])

    ax1_xz.set_xlabel('X[km]')
    ax1_xz.set_ylabel('Z[km]')

#--ax1-color-bar
    ax1_cb = divider1.append_axes('bottom', size=cb_size, pad=0.5)
    cb1 = fig.colorbar(hnd, cax=ax1_cb, orientation='horizontal')
    cb1.set_label('Stacked Local-CC Amplitude')



#--ax2: trigger grid
#--ax2_xy
    ax2_xy = fig.add_subplot(223)
    divider11 = make_axes_locatable(ax2_xy)
    hnd2 = ax2_xy.imshow(np.flipud(np.transpose(proj_grid[:,:,k_grid])),
                       extent=proj_grid.get_xy_extent(), cmap=scmap2,
                       vmin=lcc_min, vmax=lcc_max,
                       rasterized=True)
    ax2_xy.axis('tight')
    ax2_xy.set_xlim(Xmin,Xmax)
    ax2_xy.set_ylim(Ymin,Ymax)
    ax2_xy.set_xlabel('X[km]')
    ax2_xy.set_ylabel('Y[km]')
    labels = ax2_xy.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    if coord_eq:
        ax2_xy.scatter(coord_eq[0], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1, c='w')
    for sta in coord_sta:
        x_sta, y_sta = coord_sta[sta]
        ax2_xy.scatter(x_sta, y_sta, marker='^', s=sta_smbl_size, linewidths=1, c='k', alpha=0.79)
        trans = ax2_xy.transData + ax2_xy.transAxes.inverted()
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        ax2_xy.text(x_sta_ax+0.02, y_sta_ax+0.02, sta, fontsize=12, color='k', transform=ax2_xy.transAxes)
    if trigger is not None:
        t = trigger.origin_time
        if t is not None:
            t_str = '%s.%03d, ' %\
                    (t.strftime('Date: %Y-%m-%d, Time: %H:%M:%S'),
                     int(round(t.microsecond/1000.)))
        else:
            t_str = ''
        if trigger.lon:
            ax2_xy.set_title('%sLon: %.4f, Lat: %.4f, Depth: %.3f km' %
                             (t_str, trigger.lon, trigger.lat, trigger.z))
        else:
            ax2_xy.set_title('%sX: %.2f km, Y: %.2f km, Depth: %.2f km' %
                             (t_str, trigger.x, trigger.y, trigger.z))
        ax2_xy.scatter(xx_max, yy_max,
                    marker='*', s=trig_smbl_size, linewidths=1, c='g')
    ax2_xy.set_aspect('equal', 'datalim')

#--ax2_yz
    ax2_yz = divider11.append_axes('right', size=yz_size, pad=0.05, sharey=ax2_xy)
    ax2_yz.imshow(np.flipud(proj_grid[i_grid,:,:]),
                      extent=proj_grid.get_zy_extent(), cmap=scmap2,
                      vmin=lcc_min, vmax=lcc_max,
                      rasterized=True)

    if trigger is not None:
        ax2_yz.scatter(zz_max, yy_max,
                    marker='*', s=trig_smbl_size, linewidths=1, c='g')

    ax2_yz.axis('tight')
    ax2_yz.set_xlim(Zmin,Zmax)
    ax2_yz.set_ylim(Ymin,Ymax)
    ax2_yz.set_xlabel('Z[km]')
    labels = ax2_yz.get_xticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    labels = ax2_yz.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    ax2_yz.yaxis.set_visible(False)

#--ax2_xz
    ax2_xz = divider11.append_axes('bottom', size=xz_size, pad=0.05, sharex=ax2_xy)
    ax2_xz.imshow(np.flipud(np.transpose(proj_grid[:,j_grid,:])),
                     extent=proj_grid.get_xz_extent(), cmap=scmap2,
                     vmin=lcc_min, vmax=lcc_max,
                     rasterized=True)

    if trigger is not None:
        ax2_xz.scatter(xx_max, zz_max,
                    marker='*', s=trig_smbl_size, linewidths=1, c='g')

    ax2_xz.axis('tight')
    ax2_xz.set_xlim(Xmin,Xmax)
    ax2_xz.set_ylim(Zmin,Zmax)
    ax2_xz.set_ylim(ax2_xz.get_ylim()[::-1])

    ax2_xz.set_xlabel('X[km]')
    ax2_xz.set_ylabel('Z[km]')

#--ax2-color-bar
    cbx11 = divider11.append_axes('bottom', size=cb_size, pad=0.5)
    cb11=fig.colorbar(hnd2, cax=cbx11,orientation='horizontal',
                      ticks=[LTrig, LTrig+(lcc_max-LTrig)/2,lcc_max])
    cb11.set_label('Stacked Local-CC Amplitude')


#--ax3: traces
    st_plt = st.copy()
    st_plt.filter('bandpass', freqmin=fq[n22], freqmax=fq[n1],
              corners=2, zerophase=True)
    ax3 = fig.add_subplot(122)
    time += config.cut_start
    time_env += config.cut_start
    ax3.set_xlim(min(time), max(time))
    sta_y = [coord_sta[sta][1] for sta in coord_sta]
    ax3.set_ylim(min(sta_y), max(sta_y))
    ax3.set_xlabel('Time[sec]')
    ax3.set_ylabel('Y[km]')
    labels = ax3.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    trans = ax3.transData + ax3.transAxes.inverted()
    invtrans = trans.inverted()
    for sta in set(tr.stats.station for tr in st_plt):
        x_sta, y_sta = coord_sta[sta]
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        if plot_waveforms:
            # try selecting vertical component...
            try:
                tr = st_plt.select(station=sta, component='Z')[0]
            except IndexError:
                tr = None
            # otherwhise, just use the first one.
            if not tr:
                tr = st_plt.select(station=sta)[0]
            # Project signal to Axes coordinates:
            signal = tr.data/abs(tr.max())*0.05 + y_sta_ax
            xydata = np.dstack((np.zeros_like(signal), signal))[0]
            ydata = invtrans.transform(xydata)[:,1]
            ax3.plot(time, ydata, 'k', alpha=0.4, rasterized=True)
        # Project signal to Axes coordinates:
        for wave in config.wave_type:
            tr_CF = st_CF.select(station=sta, channel=wave)[0]
            signal = tr_CF.data/abs(tr_CF.max())*0.05 + y_sta_ax
            xydata = np.dstack((np.zeros_like(signal), signal))[0]
            ydata = invtrans.transform(xydata)[:,1]
            if wave == 'P':
                color = 'blue'
            if wave == 'S':
                color = 'red'
            ax3.plot(time_env, ydata, color, rasterized=True)
        ax3.text(max(time), y_sta, tr.id, fontsize=10)

        ## plotting vertical bars corresponding to LCCmax in given time window
        if trigger is not None:
            y_max = max(ydata)
            y_min = 2 * min(ydata) - y_max
            for pick in trigger.get_picks(station=sta):
                wave = pick.arrival_type
                if wave == 'P':
                    color = 'blue'
                if wave == 'S':
                    color = 'red'
                #for p_times in arrival_times[sta][wave]:
                #    LCCmax_time = p_times - st[0].stats.starttime + config.cut_start
                #    ax3.plot((LCCmax_time, LCCmax_time), (y_min, y_max), linewidth=1, color=color)
                pick_time = trigger.origin_time + pick.pick_time - st[0].stats.starttime + config.cut_start
                theor_time = trigger.origin_time + pick.theor_time - st[0].stats.starttime + config.cut_start
                if pick.pick_time != -10.:
                    ax3.plot((pick_time, pick_time), (y_min, y_max), linewidth=2.0, color=color)
                ax3.plot((theor_time, theor_time), (y_min, y_max), linewidth=2.0, color=color, linestyle='--')

    if Mtau is not None:
        for tt in Mtau:
            ax3.axvspan(t_b+config.cut_start, t_b+tt+config.cut_start,
                        facecolor='g', alpha=0.1)
    else:
        ax3.axvspan(t_b+config.cut_start, t_e+config.cut_start,
                    facecolor='g', alpha=0.1)

    note_t='CF of MBFilter; Fq= '+str(np.round(fq[n22]))+\
            '-'+str(np.round(fq[n1]))+' Hz'
    #fq_str=str(np.round(fq[n1]))+'_'+str(np.round(fq[n22]))
    ax3.set_title(note_t, fontsize=15)
    ax3.autoscale(enable=True, axis='y', tight=False)


    file_out_fig = datestr + '_t' +\
                   str('%05.1f' % (config.cut_start+t_b)) + 's_' + fq_str + '_fig.' + config.plot_format
    file_out_fig = os.path.join(out_dir, file_out_fig)
    if config.plot_format == 'pdf':
        fig.patch.set_alpha(0.0)
    # Source: http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(file_out_fig)


def plt_SummaryOut(config, grid1, st_CF, st, time_env, time, coord_sta,
                   triggers, t_bb, datestr, fq_1, fq_2,
                   coord_eq, coord_jma, file_out_fig):

    plot_waveforms = config.plot_waveforms
    ch_function = config.ch_function
    time_lag = config.time_lag

    Xmin, Xmax, Ymin, Ymax, Zmin, Zmax = grid1.get_extent()
    ratio = (Xmax - Xmin) / (Ymax - Ymin)
    fig = figure.Figure(figsize=(20, 20))
    plot_xz_size = ((Zmax - Zmin)/(Xmax - Xmin))*100
    plot_yz_size = plot_xz_size / ratio
    plot_cbar_size = 5 #percent
    xz_size = '%f %%' % plot_xz_size
    yz_size = '%f %%' % plot_yz_size
    cb_size = '%f %%' % plot_cbar_size
    sta_smbl_size = 150 / ratio
    eq_smbl_size = 200 / ratio
    trig_smbl_size = 200 / ratio

    x_trig = [ t.x for t in triggers ]
    y_trig = [ t.y for t in triggers ]
    z_trig = [ t.z for t in triggers ]

#--ax1_xy
    ax1_xy = fig.add_subplot(221)
    divider1 = make_axes_locatable(ax1_xy)
    ax1_xy.axis('tight')
    ax1_xy.set_xlim(Xmin, Xmax)
    ax1_xy.set_ylim(Ymin, Ymax)
    ax1_xy.set_xlabel('X[km]')
    ax1_xy.set_ylabel('Y[km]')
    labels = ax1_xy.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    note='Day: ' + datestr[0:6] + ',  Hour: ' + datestr[6:8]
    ax1_xy.set_title(note, fontsize=15)
    ax1_xy.scatter(x_trig, y_trig, marker='*', s=trig_smbl_size, linewidths=0.5, c='g', alpha=0.7)
    for sta in coord_sta:
        x_sta, y_sta = coord_sta[sta]
        ax1_xy.scatter(x_sta, y_sta, marker='^', s=sta_smbl_size, linewidths=1,c='k',alpha=0.79)
        trans = ax1_xy.transData + ax1_xy.transAxes.inverted()
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        ax1_xy.text(x_sta_ax+0.02, y_sta_ax+0.02, sta, fontsize=12, color='k', transform=ax1_xy.transAxes)
    if coord_eq:
        ax1_xy.scatter(coord_eq[0], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1, c='r')
    if coord_jma:
        ax1_xy.scatter(coord_jma[0], coord_jma[1], marker='o', s=eq_smbl_size, linewidths=1, c='m')
    ax1_xy.set_aspect('equal', 'datalim')

#--ax1_yz
    ax1_yz = divider1.append_axes('right', size=yz_size, pad=0.05, sharey=ax1_xy)
    ax1_yz.axis('tight')
    ax1_yz.set_xlim(Zmin,Zmax)
    ax1_yz.set_ylim(Ymin,Ymax)
    ax1_yz.set_xlabel('Z[km]')
    labels = ax1_yz.get_xticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    labels = ax1_yz.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    ax1_yz.yaxis.set_visible(False)
    ax1_yz.scatter(z_trig,y_trig,marker='*', s=trig_smbl_size, linewidths=0.5, c='g', alpha=0.7)
    if coord_eq:
        ax1_yz.scatter(coord_eq[2], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1, c='r')
    if coord_jma:
        ax1_yz.scatter(coord_jma[2], coord_jma[1], marker='o', s=eq_smbl_size, linewidths=1, c='m')

#--ax1_xz
    ax1_xz = divider1.append_axes('bottom', size=xz_size, pad=0.05, sharex=ax1_xy)
    ax1_xz.axis('tight')
    ax1_xz.set_xlim(Xmin,Xmax)
    ax1_xz.set_ylim(Zmin,Zmax)
    ax1_xz.set_ylim(ax1_xz.get_ylim()[::-1])
    ax1_xz.set_xlabel('X[km]')
    ax1_xz.set_ylabel('Z[km]')
    ax1_xz.scatter(x_trig,z_trig,marker='*', s=trig_smbl_size, linewidths=0.5, c='g', alpha=0.7)
    if coord_eq:
        ax1_xz.scatter(coord_eq[0], coord_eq[2], marker='*', s=eq_smbl_size, linewidths=1, c='r')
    if coord_jma:
        ax1_xz.scatter(coord_jma[0], coord_jma[2], marker='o', s=eq_smbl_size, linewidths=1, c='m')

#--ax1-color-bar
    # we create a color bar axis and then make it invisible,
    # just to keep the same proportions of the other plots
    ax1_cb = divider1.append_axes('bottom', size=cb_size, pad=0.5)
    ax1_cb.set_visible(False)

#--ax3: traces
    st_plt = st.copy()
    st_plt.filter('bandpass', freqmin=fq_2, freqmax=fq_1,
              corners=2, zerophase=True)
    ax3 = fig.add_subplot(122)
    time +=config.cut_start
    time_env +=config.cut_start
    ax3.set_xlim(min(time), max(time))
    sta_y = [coord_sta[sta][1] for sta in coord_sta]
    ax3.set_ylim(min(sta_y), max(sta_y))
    ax3.set_xlabel('Time[sec]')
    ax3.set_ylabel('Y[km]')
    labels = ax3.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    trans = ax3.transData + ax3.transAxes.inverted()
    invtrans = trans.inverted()
    for sta in set(tr.stats.station for tr in st_plt):
        x_sta, y_sta = coord_sta[sta]
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        if plot_waveforms:
            # try selecting vertical component...
            try:
                tr = st_plt.select(station=sta, component='Z')[0]
            except IndexError:
                tr = None
            # otherwhise, just use the first one.
            if not tr:
                tr = st_plt.select(station=sta)[0]
            # Project signal to Axes coordinates:
            signal = tr.data/abs(tr.max())*0.05 + y_sta_ax
            xydata = np.dstack((np.zeros_like(signal), signal))[0]
            ydata = invtrans.transform(xydata)[:,1]
            ax3.plot(time, ydata, 'k', alpha=0.4, rasterized=True)
        # Project signal to Axes coordinates:
        for wave in config.wave_type:
            tr_CF = st_CF.select(station=sta, channel=wave)[0]
            signal = tr_CF.data/abs(tr_CF.max())*0.05 + y_sta_ax
            xydata = np.dstack((np.zeros_like(signal), signal))[0]
            ydata = invtrans.transform(xydata)[:,1]
            if wave == 'P':
                color = 'blue'
            if wave == 'S':
                color = 'red'
            ax3.plot(time_env, ydata, color, rasterized=True)
        ax3.text(max(time), y_sta, tr.id, fontsize=10)

    note=ch_function+' of MBFilter; Fq. range: '+str(np.round(fq_1))+\
            '-'+str(np.round(fq_2))+' Hz'
    ax3.set_title(note, fontsize=15)
    ##[ax3.axvline(center_trigWin[m],linewidth=2, color='r',alpha=0.2)\
    ##     for m in xrange(len(beg_trigWin))]

    for t in triggers:
        ax3.axvspan(t.beg_win, t.end_win,
                    facecolor='r', alpha=0.1)

    ax3.axvline(t_bb[0]+config.cut_start,linewidth=1, color='b',alpha=0.9)
    ax3.axvline(t_bb[-1]+time_lag+config.cut_start,linewidth=1, color='b',alpha=0.9)
    ax3.autoscale(enable=True, axis='y', tight=False)

    if config.plot_format == 'pdf':
        fig.patch.set_alpha(0.0)
    # Source: http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(file_out_fig)
