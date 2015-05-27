import os
import numpy as np
import pylab
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.patheffects as path_effects


def bp_plot(config, proj_grid,
            coord_eq, t_begin, t_end,
            coord_sta,
            st, st_CF,
            freqs, trigger,
            arrival_times=None,
            unused_CF=None,
            Mtau=None,
            async_plotter=None):

    if trigger is not None:
        LTrig = trigger.trigger_level
        lcc_max = proj_grid.max()
    else:
        LTrig = config.trigger
        lcc_max = 1
    lcc_min = LTrig
    out_dir = config.out_dir
    plot_waveforms = config.plot_waveforms

    Xmin, Xmax, Ymin, Ymax, Zmin, Zmax = proj_grid.get_extent()
    ratio = (Xmax - Xmin) / (Ymax - Ymin)
    sta_smbl_size = 150 / ratio
    eq_smbl_size = 200 / ratio
    trig_smbl_size = 200 / ratio

    scmap = config.scmap
    scmap2 = pylab.cm.jet
    scmap2.set_under('w', LTrig)

    if trigger is not None:
        max_ijk = map(int, (trigger.i, trigger.j, trigger.k))
        xx_max, yy_max, zz_max = trigger.x, trigger.y, trigger.z
    elif coord_eq:
        x_eq, y_eq, z_eq = coord_eq
        max_ijk = proj_grid.get_ijk(x_eq[0], y_eq[0], z_eq[0])
    else:
        max_ijk = proj_grid.get_ijk_max()

    try:
        box_idx = proj_grid.box_idx
        i1, i2, j1, j2, k1, k2 = box_idx
    except AttributeError:
        box_idx = None

#--figure
    fig = figure.Figure(figsize=(20, 20))

#--ax1: stacked grid
    ax1_xy = fig.add_subplot(221)
    axes, cb1 = proj_grid.plot(max_ijk, handle=True, ax_xy=ax1_xy, vmin=0, cmap=scmap)
    cb1.set_label('Stacked Local-CC Amplitude')

    ax1_xy, ax1_xz, ax1_yz = axes
    ax1_xy.set_xlabel('X[km]')
    ax1_xy.set_ylabel('Y[km]')
    ax1_yz.set_xlabel('Z[km]')
    ax1_xz.set_xlabel('X[km]')
    ax1_xz.set_ylabel('Z[km]')

    tt1 = st[0].stats.starttime + t_begin
    tt2 = st[0].stats.starttime + t_end
    ax1_xy.set_title('Date: %s, Time: %s.%03d - %s.%03d (%s - %s s)' %
                     (tt1.date,
                      tt1.strftime('%H:%M:%S'),
                      int(round(tt1.microsecond/1000.)),
                      tt2.strftime('%H:%M:%S'),
                      int(round(tt2.microsecond/1000.)),
                      t_begin + config.cut_start,
                      t_end + config.cut_start))
    if coord_eq:
        ax1_xy.scatter(coord_eq[0], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1,c='w')
    for sta in coord_sta:
        x_sta, y_sta = coord_sta[sta]
        ax1_xy.scatter(x_sta, y_sta, marker='^', s=sta_smbl_size, linewidths=1, c='w', alpha=0.79)
        trans = ax1_xy.transData + ax1_xy.transAxes.inverted()
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        ax1_xy.text(x_sta_ax+0.02, y_sta_ax+0.02, sta, fontsize=12, color='w',
                    transform=ax1_xy.transAxes,
                    path_effects=[path_effects.withStroke(linewidth=2, foreground='k')])

    if trigger is not None:
        ax1_xy.scatter(xx_max, yy_max,
                       marker='*', s=trig_smbl_size, linewidths=1, c='g')
        ax1_yz.scatter(zz_max, yy_max,
                       marker='*', s=trig_smbl_size, linewidths=1, c='g')
        ax1_xz.scatter(xx_max, zz_max,
                       marker='*', s=trig_smbl_size, linewidths=1, c='g')
    if proj_grid.ellipsoid is not None:
        proj_grid.plot_ellipsoid(axes, proj_grid.ellipsoid, proj_grid.xyz_mean)
    if box_idx is not None:
        x_xy = np.array((i1, i2, i2, i1, i1)) * proj_grid.dx + proj_grid.x_orig
        y_xy = np.array((j1, j1, j2, j2, j1)) * proj_grid.dy + proj_grid.y_orig
        z_yz = np.array((k1, k2, k2, k1, k1)) * proj_grid.dz + proj_grid.z_orig
        y_yz = np.array((j1, j1, j2, j2, j1)) * proj_grid.dy + proj_grid.y_orig
        x_xz = np.array((i1, i2, i2, i1, i1)) * proj_grid.dx + proj_grid.x_orig
        z_xz = np.array((k1, k1, k2, k2, k1)) * proj_grid.dz + proj_grid.z_orig
        ax1_xy.plot(x_xy, y_xy, color='w')
        ax1_yz.plot(z_yz, y_yz, color='w')
        ax1_xz.plot(x_xz, z_xz, color='w')


#--ax2: trigger grid
    if trigger is not None:
        ax2_xy = fig.add_subplot(223)
        axes, cb11 = proj_grid.plot(max_ijk, handle=True, ax_xy=ax2_xy, vmin=lcc_min, vmax=lcc_max, cmap=scmap2)
        cb11.set_label('Stacked Local-CC Amplitude')
        cb11.set_ticks([LTrig, LTrig+(lcc_max-LTrig)/2, lcc_max])

        ax2_xy, ax2_xz, ax2_yz = axes
        ax2_xy.set_xlabel('X[km]')
        ax2_xy.set_ylabel('Y[km]')
        ax2_yz.set_xlabel('Z[km]')
        ax2_xz.set_xlabel('X[km]')
        ax2_xz.set_ylabel('Z[km]')

        if coord_eq:
            ax2_xy.scatter(coord_eq[0], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1, c='w')
        for sta in coord_sta:
            x_sta, y_sta = coord_sta[sta]
            ax2_xy.scatter(x_sta, y_sta, marker='^', s=sta_smbl_size, linewidths=1, c='k', alpha=0.79)
            trans = ax2_xy.transData + ax2_xy.transAxes.inverted()
            x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
            ax2_xy.text(x_sta_ax+0.02, y_sta_ax+0.02, sta, fontsize=12, color='k', transform=ax2_xy.transAxes)

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
        ax2_yz.scatter(zz_max, yy_max,
                       marker='*', s=trig_smbl_size, linewidths=1, c='g')
        ax2_xz.scatter(xx_max, zz_max,
                       marker='*', s=trig_smbl_size, linewidths=1, c='g')


#--ax3: traces
    st_plt = st.copy()
    st_plt.detrend(type='constant')
    st_plt.filter('bandpass', freqmin=freqs[-1], freqmax=freqs[0],
                  corners=2, zerophase=True)
    st_CF_plt = st_CF.copy()
    if config.plot_time_win_size is not None:
        time_win_size = config.plot_time_win_size
        t_extra = (time_win_size - (t_end - t_begin)) / 2
        st_plt.trim(config.starttime + t_begin - t_extra,
                    config.starttime + t_end + t_extra)
        st_CF_plt.trim(config.starttime + t_begin - t_extra,
                       config.starttime + t_end + t_extra)

    ax3 = fig.add_subplot(122)
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
            time = np.arange(tr.stats.npts) / tr.stats.sampling_rate
            if config.plot_time_win_size is not None:
                time += config.cut_start + t_begin - t_extra
                # remove negative times by shifting
                if time.min() < 0:
                    time -= time.min()
            else:
                time += config.cut_start
            ax3.plot(time, ydata, 'k', alpha=0.4, rasterized=True)
        # Project signal to Axes coordinates:
        for wave in config.wave_type:
            tr_CF = st_CF_plt.select(station=sta, channel=wave)[0]
            signal = tr_CF.data/abs(tr_CF.max())*0.05 + y_sta_ax
            xydata = np.dstack((np.zeros_like(signal), signal))[0]
            ydata = invtrans.transform(xydata)[:,1]
            if wave == 'P':
                color = 'blue'
            if wave == 'S':
                color = 'red'
            time_CF = np.arange(tr_CF.stats.npts) / tr_CF.stats.sampling_rate
            if config.plot_time_win_size is not None:
                time_CF += config.cut_start + t_begin - t_extra
                # remove negative times by shifting
                if time_CF.min() < 0:
                    time_CF -= time_CF.min()
            else:
                time_CF += config.cut_start
            if '%s.%s' % (sta, wave) in unused_CF:
                linestyle = '--'
            else:
                linestyle = '-'
            ax3.plot(time_CF, ydata,
                     color=color, linestyle=linestyle,
                     rasterized=True)
        ax3.set_xlim(min(time), max(time))
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
                theor_time = trigger.origin_time + pick.theor_time - st[0].stats.starttime + config.cut_start
                if pick.pick_time != -10.:
                    pick_time = trigger.origin_time + pick.pick_time - st[0].stats.starttime + config.cut_start
                    ax3.plot((pick_time, pick_time), (y_min, y_max), linewidth=2.0, color=color)
                ax3.plot((theor_time, theor_time), (y_min, y_max), linewidth=2.0, color=color, linestyle='--')

    if Mtau is not None:
        for tt in Mtau:
            ax3.axvspan(t_begin+config.cut_start, t_begin+tt+config.cut_start,
                        facecolor='g', alpha=0.1)
    else:
        ax3.axvspan(t_begin+config.cut_start, t_end+config.cut_start,
                    facecolor='g', alpha=0.1)

    note_t = 'CF of MBFilter; Fq= ' + str(np.round(freqs[-1])) +\
             '-' + str(np.round(freqs[0])) + ' Hz'
    ax3.set_title(note_t, fontsize=15)
    ax3.autoscale(enable=True, axis='y', tight=False)


    fq_str = str(np.round(freqs[-1])) + '_' + str(np.round(freqs[0]))
    datestr = st[0].stats.starttime.strftime('%y%m%d%H')

    file_out_fig = datestr + '_t' +\
                   str('%06.1f' % (config.cut_start+t_begin)) + 's_' + fq_str + '_fig.' + config.plot_format
    file_out_fig = os.path.join(out_dir, file_out_fig)
    if config.plot_format == 'pdf':
        fig.patch.set_alpha(0.0)
    # Source: http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
    canvas = FigureCanvasAgg(fig)
    if async_plotter is not None:
        async_plotter.save(canvas, file_out_fig)
    else:
        canvas.print_figure(file_out_fig)
    fig.clf()


def plt_SummaryOut(config, grid1, st_CF, st, coord_sta,
                   triggers, t_bb, datestr, fq_1, fq_2,
                   coord_eq, coord_jma, file_out_fig):

    plot_waveforms = config.plot_waveforms
    ch_function = config.ch_function
    time_lag = config.time_lag

    Xmin, Xmax, Ymin, Ymax, Zmin, Zmax = grid1.get_extent()
    ratio = (Xmax - Xmin) / (Ymax - Ymin)
    sta_smbl_size = 150 / ratio
    eq_smbl_size = 200 / ratio
    trig_smbl_size = 200 / ratio

    x_trig = [ t.x for t in triggers ]
    y_trig = [ t.y for t in triggers ]
    z_trig = [ t.z for t in triggers ]

#-- figure
    fig = figure.Figure(figsize=(20, 20))

#--ax1: summary plot
    ax1_xy = fig.add_subplot(221)
    ax1_xy, ax1_xz, ax1_yz, ax1_cb = grid1.get_plot_axes(figure, ax1_xy)
    ax1_cb.set_visible(False)

    ax1_xy.set_xlabel('X[km]')
    ax1_xy.set_ylabel('Y[km]')
    ax1_yz.set_xlabel('Z[km]')
    ax1_xz.set_xlabel('X[km]')
    ax1_xz.set_ylabel('Z[km]')

    note = 'Day: ' + datestr[0:6] + ',  Hour: ' + datestr[6:8]
    ax1_xy.set_title(note, fontsize=15)
    ax1_xy.scatter(x_trig, y_trig, marker='*', s=trig_smbl_size, linewidths=0.5, c='g', alpha=0.7)
    ax1_yz.scatter(z_trig, y_trig, marker='*', s=trig_smbl_size, linewidths=0.5, c='g', alpha=0.7)
    ax1_xz.scatter(x_trig, z_trig, marker='*', s=trig_smbl_size, linewidths=0.5, c='g', alpha=0.7)
    for sta in coord_sta:
        x_sta, y_sta = coord_sta[sta]
        ax1_xy.scatter(x_sta, y_sta, marker='^', s=sta_smbl_size, linewidths=1,c='k',alpha=0.79)
        trans = ax1_xy.transData + ax1_xy.transAxes.inverted()
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        ax1_xy.text(x_sta_ax+0.02, y_sta_ax+0.02, sta, fontsize=12, color='k', transform=ax1_xy.transAxes)
    if coord_eq:
        ax1_xy.scatter(coord_eq[0], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1, c='r')
        ax1_yz.scatter(coord_eq[2], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1, c='r')
        ax1_xz.scatter(coord_eq[0], coord_eq[2], marker='*', s=eq_smbl_size, linewidths=1, c='r')
    if coord_jma:
        ax1_xy.scatter(coord_jma[0], coord_jma[1], marker='o', s=eq_smbl_size, linewidths=1, c='m')
        ax1_yz.scatter(coord_jma[2], coord_jma[1], marker='o', s=eq_smbl_size, linewidths=1, c='m')
        ax1_xz.scatter(coord_jma[0], coord_jma[2], marker='o', s=eq_smbl_size, linewidths=1, c='m')

#--ax3: traces
    st_plt = st.copy()
    st_plt.detrend(type='constant')
    st_plt.filter('bandpass', freqmin=fq_2, freqmax=fq_1,
                  corners=2, zerophase=True)
    ax3 = fig.add_subplot(122)
    time = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate
    time += config.cut_start
    time_env = np.arange(st_CF[0].stats.npts) / st_CF[0].stats.sampling_rate
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

    note = ch_function + ' of MBFilter; Fq. range: ' + str(np.round(fq_1)) +\
           '-' + str(np.round(fq_2)) + ' Hz'
    ax3.set_title(note, fontsize=15)

    for t in triggers:
        ax3.axvspan(t.beg_win, t.end_win,
                    facecolor='r', alpha=0.1)

    ax3.axvline(t_bb[0]+config.cut_start, linewidth=1, color='b', alpha=0.9)
    ax3.axvline(t_bb[-1]+time_lag+config.cut_start, linewidth=1, color='b', alpha=0.9)
    ax3.autoscale(enable=True, axis='y', tight=False)

    if config.plot_format == 'pdf':
        fig.patch.set_alpha(0.0)
    # Source: http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(file_out_fig)
    fig.clf()
