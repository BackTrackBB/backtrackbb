import os
import numpy as np
import pylab
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.axes_grid1 import make_axes_locatable

def bp_plot(config, grid1, proj_grid, comb_sta,
        coord_eq, t_b, t_e, datestr, fq_str,
        extent_grd, extent_yz, extent_xz,
        coord_sta, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
        st, sta, st_CF,
        time, time_env,
        fq, n1, n22):

    LTrig = config.trigger
    lcc_max = config.lcc_max
    time_lag = config.time_lag
    out_dir = config.out_dir
    scmap = config.scmap
    plot_waveforms = config.plot_waveforms

    ratio = (Xmax - Xmin) / (Ymax - Ymin)
    fig_size_x = 18
    fig_size_y = fig_size_x / ratio
    fig = figure.Figure(figsize=(fig_size_x, fig_size_y))
    plot_xz_size = 25 #percent
    plot_yz_size = plot_xz_size / ratio
    xz_size = '%f %%' % plot_xz_size
    yz_size = '%f %%' % plot_yz_size
    sta_smbl_size = 250 / ratio
    eq_smbl_size = 300 / ratio
    trig_smbl_size = 200 / ratio

    #grid_max= proj_grid/len(comb_sta)
    grid_max= proj_grid/1

    scmap2 = pylab.cm.jet
    scmap2.set_under('w', LTrig)
    lcc_min = LTrig

    Max_grid_max =np.where(grid_max == np.max(grid_max))
    x_max = Max_grid_max[0][0]
    y_max = Max_grid_max[1][0]
    z_max = Max_grid_max[2][0]

    extent_zy = (extent_yz[2],extent_yz[3],extent_yz[0],extent_yz[1])

#--ax1: stacked grid
#--ax1_xy
    ax1_xy = fig.add_subplot(221)
    divider1 = make_axes_locatable(ax1_xy)

    if grid_max[x_max,y_max,z_max] >= LTrig:
        z_grid = z_max
        y_grid = y_max
        x_grid = x_max
    elif coord_eq:
        x_eq, y_eq, z_eq = coord_eq
        n_y = np.argwhere(grid1.y_array >= y_eq[0])[0][0]
        n_x = np.argwhere(grid1.x_array >= x_eq[0])[0][0]
        n_z = np.argwhere(grid1.z_array >= z_eq)[0][0]
        z_grid = n_z
        y_grid = n_y
        x_grid = n_x
    else:
        x_grid = int(proj_grid.shape[0] / 2)
        y_grid = int(proj_grid.shape[1] / 2)
        z_grid = int(proj_grid.shape[2] / 2)
    hnd=ax1_xy.imshow(np.flipud(np.transpose(grid_max[:,:,z_grid])),
                             extent=extent_grd,cmap=scmap,rasterized=True)
    ax1_xy.axis('tight')
    ax1_xy.set_xlim(Xmin,Xmax)
    ax1_xy.set_ylim(Ymin,Ymax)
    ax1_xy.set_xlabel('X[km]')
    ax1_xy.set_ylabel('Y[km]')
    labels = ax1_xy.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    tt1 = st[0].stats.starttime+t_b
    tt2 = st[0].stats.starttime+t_e
    ax1_xy.set_title('Date: '+str(tt1.date)+', Time: '+str(tt1.time)+\
                  ' - '+str(tt2.time)+'  ('+str(t_b)+' - '+str(t_e)+'  sec),'+
                  ' depth='+str(int(grid1.z_array[z_grid]))+' [km]')
    if coord_eq:
        ax1_xy.scatter(coord_eq[0], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1,c='w')
    for sta in coord_sta:
        x_sta, y_sta = coord_sta[sta]
        ax1_xy.scatter(x_sta, y_sta, marker='^', s=sta_smbl_size, linewidths=1, c='k', alpha=0.79)
        trans = ax1_xy.transData + ax1_xy.transAxes.inverted()
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        ax1_xy.text(x_sta_ax+0.02, y_sta_ax+0.02, sta, fontsize=12, color='k', transform=ax1_xy.transAxes)

    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax1_xy.scatter(grid1.x_array[x_max],grid1.y_array[y_max],
                    marker='*', s=trig_smbl_size, linewidths=1,c='g')

#--ax1_yz
    ax1_yz = divider1.append_axes('right', size=yz_size, pad=0.05, sharey=ax1_xy)
    ax1_yz.imshow(np.flipud(grid_max[x_grid,:,:]),
                     extent=extent_zy,cmap=scmap,rasterized=True)

    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax1_yz.scatter(grid1.z_array[z_max],grid1.y_array[y_max],
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
    ax1_xz.imshow(np.flipud(np.transpose(grid_max[:,y_grid,:])),
                     extent=extent_xz,cmap=scmap,rasterized=True)

    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax1_xz.scatter(grid1.x_array[x_max],grid1.z_array[z_max],
                    marker='*', s=trig_smbl_size, linewidths=1, c='g')                      

    ax1_xz.axis('tight')
    ax1_xz.set_xlim(Xmin,Xmax)
    ax1_xz.set_ylim(Zmin,Zmax)
    ax1_xz.set_ylim(ax1_xz.get_ylim()[::-1])

    ax1_xz.set_xlabel('X[km]')
    ax1_xz.set_ylabel('Z[km]')

#--ax1-color-bar
    ax1_cb = divider1.append_axes('bottom', size='5%', pad=0.5)
    cb1 = fig.colorbar(hnd, cax=ax1_cb, orientation='horizontal')
    cb1.set_label('Stacked Local-CC Amplitude')



#--ax2: trigger grid
#--ax2_xy
    ax2_xy = fig.add_subplot(223)
    divider11 = make_axes_locatable(ax2_xy)
    hnd2 = ax2_xy.imshow(np.flipud(np.transpose(grid_max[:,:,z_grid])),
                       extent=extent_grd,cmap=scmap2,
                       vmin=lcc_min,vmax=lcc_max,
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
    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax2_xy.scatter(grid1.x_array[x_max],grid1.y_array[y_max],
                    marker='*', s=trig_smbl_size, linewidths=1, c='g')

#--ax2_yz
    ax2_yz = divider11.append_axes('right', size=yz_size, pad=0.05, sharey=ax2_xy)
    ax2_yz.imshow(np.flipud(grid_max[x_grid,:,:]),
                      extent=extent_zy,cmap=scmap2,vmin=lcc_min,vmax=lcc_max,
                      rasterized=True)

    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax2_yz.scatter(grid1.z_array[z_max],grid1.y_array[y_max],
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
    ax2_xz.imshow(np.flipud(np.transpose(grid_max[:,y_grid,:])),
                     extent=extent_xz,cmap=scmap2,
                      vmin=lcc_min,vmax=lcc_max,
                      rasterized=True)

    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax2_xz.scatter(grid1.x_array[x_max],grid1.z_array[z_max],
                    marker='*', s=trig_smbl_size, linewidths=1, c='g')

    ax2_xz.axis('tight')
    ax2_xz.set_xlim(Xmin,Xmax)
    ax2_xz.set_ylim(Zmin,Zmax)
    ax2_xz.set_ylim(ax2_xz.get_ylim()[::-1])

    ax2_xz.set_xlabel('X[km]')
    ax2_xz.set_ylabel('Z[km]')

#--ax2-color-bar
    cbx11 = divider11.append_axes('bottom', size='5%', pad=0.5)
    cb11=fig.colorbar(hnd2, cax=cbx11,orientation='horizontal',
                      ticks=[LTrig, LTrig+(lcc_max-LTrig)/2,lcc_max])
    cb11.set_label('Stacked Local-CC Amplitude')


#---##    
    #ax5 = fig.add_subplot(gs[6:10, 5:9])
    #ax5.axis('tight')
    #ax5.set_title('Areas of Stacked Normalized Local-CC >'+str(LTrig))
    #ax5.set_xlim(Xmin,Xmax)
    #ax5.set_ylim(Ymin,Ymax)
    #ax5.set_xlabel('X[km]')
    #ax5.set_ylabel('Y[km]')
    #labels = ax5.get_yticklabels()
    #pylab.setp(labels, rotation=90, fontsize=12)
    #if len(grid_max[grid_max>=LTrig])>1:  
    #    ax5.scatter(grid1.x_array[x_max],grid1.y_array[y_max],
    #                marker='*', s = 200, linewidths=1,c='g')
    #for sta in coord_sta:
    #    x_sta, y_sta = coord_sta[sta]
    #    ax5.scatter(x_sta, y_sta, marker='^', s=250, linewidths=1, c='k', alpha=0.79)
    #    trans = ax5.transData + ax5.transAxes.inverted()
    #    x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
    #    ax5.text(x_sta_ax+0.02, y_sta_ax+0.02, sta, fontsize=12, color='k', transform=ax5.transAxes)
    #if coord_eq:
    #    ax5.scatter(coord_eq[0], coord_eq[1], marker='*', s = 300, linewidths=1,c='w')


#--ax3: traces
    ax3 = fig.add_subplot(122)
    ax3.set_xlim(0, max(time))
    sta_y = [coord_sta[sta][1] for sta in coord_sta]
    ax3.set_ylim(min(sta_y), max(sta_y))
    ax3.set_xlabel('Time[sec]')
    ax3.set_ylabel('Y[km]')
    labels = ax3.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    trans = ax3.transData + ax3.transAxes.inverted()
    invtrans = trans.inverted()
    for Record, CH_fct in zip(st, st_CF):
        sta = Record.stats.station
        x_sta, y_sta = coord_sta[sta]
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        if plot_waveforms:
            # Project signal to Axes coordinates:
            signal = Record.data/abs(Record.max())*0.05 + y_sta_ax
            xydata = np.dstack((np.zeros_like(signal), signal))[0]
            ydata = invtrans.transform(xydata)[:,1]
            ax3.plot(time, ydata, 'k', alpha=0.4, rasterized=True)
        # Project signal to Axes coordinates:
        signal = CH_fct.data/abs(CH_fct.max())*0.05 + y_sta_ax
        xydata = np.dstack((np.zeros_like(signal), signal))[0]
        ydata = invtrans.transform(xydata)[:,1]
        ax3.plot(time_env, ydata, 'k', rasterized=True)
        ax3.text(max(time), y_sta, Record.id, fontsize=10)
    ax3.axvspan(t_b, t_e, facecolor='g', alpha=0.2)
    note_t='CF of MBFilter; Fq= '+str(np.round(fq[n22]))+\
            '-'+str(np.round(fq[n1]))+' Hz'
    #fq_str=str(np.round(fq[n1]))+'_'+str(np.round(fq[n22]))
    ax3.set_title(note_t, fontsize=15)
    if len(grid_max[grid_max >= LTrig]) > 1:
        ax3.axvline(t_b+time_lag/2,linewidth=4, color='r',alpha=0.15)
    ax3.autoscale(enable=True, axis='y', tight=False)

    
    file_out_fig = datestr + '_t' +\
                   str('%05.1f' % (t_b)) + 's_' + fq_str + '_fig.' + config.plot_format
    file_out_fig = os.path.join(out_dir, file_out_fig)
    if config.plot_format == 'pdf':
        fig.patch.set_alpha(0.0)
    # Source: http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(file_out_fig)
    
    
def plt_SummaryOut(config, st_CF, st, time_env, time, coord_sta,
                   x_trig, y_trig, z_trig, beg_trigWin, end_trigWin, center_trigWin, t_bb,
                   Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, datestr,
                   fq_1, fq_2,
                   coord_eq, coord_jma, file_out_fig):

    plot_waveforms = config.plot_waveforms
    ch_function = config.ch_function
    time_lag = config.time_lag
    
    ratio = (Xmax - Xmin) / (Ymax - Ymin)
    fig_size_x = 18
    fig_size_y = fig_size_x / ratio
    fig = figure.Figure(figsize=(fig_size_x, fig_size_y))
    plot_xz_size = 25 #percent
    plot_yz_size = plot_xz_size / ratio
    xz_size = '%f %%' % plot_xz_size
    yz_size = '%f %%' % plot_yz_size
    sta_smbl_size = 250 / ratio
    eq_smbl_size = 300 / ratio
    trig_smbl_size = 80 / ratio


#--ax1_xy
    ax1_xy = fig.add_subplot(221)
    divider1 = make_axes_locatable(ax1_xy)
    ax1_xy.axis('tight')
    ax1_xy.set_xlim(Xmin,Xmax)
    ax1_xy.set_ylim(Ymin,Ymax)
    ax1_xy.set_xlabel('X[km]')
    ax1_xy.set_ylabel('Y[km]')
    labels = ax1_xy.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    note='Day: ' + datestr[0:6] + ',  Hour: ' + datestr[6:8]
    ax1_xy.text(Xmin,Ymax,note,fontsize=15)
    ax1_xy.scatter(x_trig, y_trig,marker='*', s=trig_smbl_size, linewidths=0.5, c='g', alpha=0.7)    
    for sta in coord_sta:
        x_sta, y_sta = coord_sta[sta]
        ax1_xy.scatter(x_sta, y_sta,marker='^', s=sta_smbl_size, linewidths=1,c='k',alpha=0.79)
        trans = ax1_xy.transData + ax1_xy.transAxes.inverted()
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        ax1_xy.text(x_sta_ax+0.02, y_sta_ax+0.02, sta, fontsize=12, color='k', transform=ax1_xy.transAxes)
    if coord_eq:
        ax1_xy.scatter(coord_eq[0], coord_eq[1], marker='*', s=eq_smbl_size, linewidths=1, c='r')
    if coord_jma:
        ax1_xy.scatter(coord_jma[0], coord_jma[1], marker='o', s=eq_smbl_size, linewidths=1, c='m')

#--ax1_yz
    ax1_yz = divider1.append_axes('right', size=yz_size, pad=0.05, sharey=ax1_xy)
    ax1_yz.axis('tight')
    ax1_yz.set_xlim(Zmin,Zmax)
    ax1_yz.set_ylim(Ymin,Ymax)
    ax1_yz.set_xlabel('Z[km]')
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
    labels = ax1_xz.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    ax1_xz.scatter(x_trig,z_trig,marker='*', s=trig_smbl_size, linewidths=0.5, c='g', alpha=0.7)
    if coord_eq:
        ax1_xz.scatter(coord_eq[0], coord_eq[2], marker='*', s=eq_smbl_size, linewidths=1, c='r')
    if coord_jma:
        ax1_xz.scatter(coord_jma[0], coord_jma[2], marker='o', s=eq_smbl_size, linewidths=1, c='m')

#--ax3: traces
    ax3 = fig.add_subplot(122)
    ax3.set_xlim(0, max(time))
    sta_y = [coord_sta[sta][1] for sta in coord_sta]
    ax3.set_ylim(min(sta_y), max(sta_y))
    ax3.set_xlabel('Time[sec]')
    ax3.set_ylabel('Y[km]')
    labels = ax3.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    trans = ax3.transData + ax3.transAxes.inverted()
    invtrans = trans.inverted()
    for Record, CH_fct in zip(st, st_CF):
        sta = Record.stats.station
        x_sta, y_sta = coord_sta[sta]
        x_sta_ax, y_sta_ax = trans.transform((x_sta, y_sta))
        if plot_waveforms:
            # Project signal to Axes coordinates:
            signal = Record.data/abs(Record.max())*0.05 + y_sta_ax
            xydata = np.dstack((np.zeros_like(signal), signal))[0]
            ydata = invtrans.transform(xydata)[:,1]
            ax3.plot(time, ydata, 'k', alpha=0.4, rasterized=True)
        # Project signal to Axes coordinates:
        signal = CH_fct.data/abs(CH_fct.max())*0.05 + y_sta_ax
        xydata = np.dstack((np.zeros_like(signal), signal))[0]
        ydata = invtrans.transform(xydata)[:,1]
        ax3.plot(time_env, ydata, 'k', rasterized=True)
        ax3.text(max(time), y_sta, Record.id, fontsize=10)

    note=ch_function+' of MBFilter; Fq. range: '+str(np.round(fq_1))+\
            '-'+str(np.round(fq_2))+' Hz'
    ax3.set_title(note, fontsize=15)
    ##[ax3.axvline(center_trigWin[m],linewidth=2, color='r',alpha=0.2)\
    ##     for m in xrange(len(beg_trigWin))]

    [ax3.axvspan(beg_trigWin[m],end_trigWin[m],facecolor='r',alpha=0.1)\
         for m in xrange(len(beg_trigWin))]

    ax3.axvline(t_bb[0],linewidth=1, color='b',alpha=0.9)
    ax3.axvline(t_bb[-1]+time_lag,linewidth=1, color='b',alpha=0.9)
    ax3.autoscale(enable=True, axis='y', tight=False)

    if config.plot_format == 'pdf':
        fig.patch.set_alpha(0.0)
    # Source: http://www.dalkescientific.com/writings/diary/archive/2005/04/23/matplotlib_without_gui.html
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(file_out_fig)


def bp_plot_pdf(grid1, grid_max, proj_pdf, x_eq, y_eq,z_eq, LTrig,
        t_b, t_e, out_dir, datestr, fq_str,
        extent_grd, extent_yz, extent_xz,
        x_sta, y_sta,
        Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
        st, scmap, lcc_min, lcc_max,
        sta, st_CF,
        time, time_env, time_lag,
        plot_waveforms, fq,
        n1, n22):
    
    scmap = pylab.cm.jet
    scmap.set_under('w', LTrig)
    lcc_min = LTrig

    fig = figure.Figure(figsize=(18,17))
    #grid_max= proj_grid/len(comb_sta)
    #grid_max= proj_grid/1

    Max_grid_max =np.where(grid_max == np.max(grid_max))
    x_max = Max_grid_max[0][0]
    y_max = Max_grid_max[1][0]
    z_max = Max_grid_max[2][0]

    n_y = np.argwhere(grid1.y_array >= y_eq[0])[0][0]
    n_x = np.argwhere(grid1.x_array >= x_eq[0])[0][0]
    n_z = np.argwhere(grid1.z_array >= z_eq)[0][0]

#--plotting horizontal projection of the stacked grid----------------------
    ax1 = fig.add_axes([0.03,0.635,0.4, 0.345])
    cbx1 = fig.add_axes([0.03, 0.515, 0.4, 0.01])
    if grid_max[x_max,y_max,z_max] >= LTrig:
        z_grid=z_max
        y_grid=y_max
        x_grid=x_max
        print str(grid1.x_array[x_max])+'  '+str(grid1.y_array[y_max])+\
              ' '+str(grid1.z_array[z_max])+'  '+str(t_b)+' '+str(t_e)       
    else:
        z_grid=n_z
        y_grid=n_y
        x_grid=n_x
        
    hnd=ax1.imshow(np.flipud(np.transpose(grid_max[:,:,z_grid])),
                   extent=extent_grd,cmap=scmap,
                   vmin=lcc_min,vmax=lcc_max,
                   rasterized=True)
    cb1=ax1.colorbar(hnd, cax=cbx1,orientation='horizontal')
    ax1.axis('tight')
    ax1.scatter(x_eq,y_eq,marker='*', s = 300, linewidths=1,c='w')
    ax1.scatter(x_sta,y_sta,marker='^', s = 250, linewidths=1,c='k',alpha=0.79)

    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax1.scatter(grid1.x_array[x_max],grid1.y_array[y_max],
                    marker='*', s = 200, linewidths=1,c='g')
    ax1.set_xlim(Xmin,Xmax)
    ax1.set_ylim(Ymin,Ymax)
    ax1.set_xlabel('X[km]')
    ax1.set_ylabel('Y[km]')
    cb1.set_label('Stacked Local-CC Amplitude')
    labels = ax1.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    tt1 = st[0].stats.starttime+t_b
    tt2 = st[0].stats.starttime+t_e
    
    ax1.set_title('Date: '+str(tt1.date)+', Time: '+str(tt1.time)+\
                  ' - '+str(tt2.time)+'  ('+str(t_b)+' - '+str(t_e)+'  sec),'+
                  ' depth='+str(int(grid1.z_array[z_grid]))+' [km]')
    

    ax11 = fig.add_axes([0.03,0.136,0.4, 0.345])
    cbx11 = fig.add_axes([0.03, 0.0155, 0.4, 0.01])
    hnd2 = ax11.imshow(np.flipud(np.transpose(proj_pdf[:,:,z_grid])),
                       extent=extent_grd,cmap=scmap,
                       vmin=lcc_min,vmax=lcc_max,
                       rasterized=True)
    ax11.axis('tight')
    cb11=ax11.colorbar(hnd2, cax=cbx11,orientation='horizontal')
    ax11.scatter(x_eq,y_eq,marker='*', s = 300, linewidths=1,c='w')
    ax11.scatter(x_sta,y_sta,marker='^', s = 250, linewidths=1,c='k',alpha=0.79)

    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax11.scatter(grid1.x_array[x_max],grid1.y_array[y_max],
                    marker='*', s = 200, linewidths=1,c='g')
    
    ax11.set_xlim(Xmin,Xmax)
    ax11.set_ylim(Ymin,Ymax)
    ax11.set_xlabel('X[km]')
    ax11.set_ylabel('Y[km]')
    cb11.set_label('Stacked Local-CC Amplitude')
    labels = ax11.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)

#---    
    ax5 = fig.add_axes([0.545,0.04,0.4, 0.345])
    if len(grid_max[grid_max>=LTrig])>1:  
        ax5.scatter(grid1.x_array[x_max],grid1.y_array[y_max],
                    marker='*', s = 200, linewidths=1,c='g')
    ax5.scatter(x_sta,y_sta,marker='^', s = 250, linewidths=1,c='k',alpha=0.79)
    ax5.scatter(x_eq,y_eq,marker='*', s = 300, linewidths=1,c='w')
    ax5.set_title('Areas of Stacked Normalized Local-CC >'+str(LTrig))
    ax5.axis('tight')
    ax5.set_xlim(Xmin,Xmax)
    ax5.set_ylim(Ymin,Ymax)
    ax5.set_xlabel('X[km]')
    ax5.set_ylabel('Y[km]')
    labels = ax5.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)

    
#--plotting vertical projections of the stacked grid
    ax2 = fig.add_axes([0.45,0.63,0.05,0.345])
    ax22 = fig.add_axes([0.45,0.136,0.05,0.345])

    extent_zy = (extent_yz[2],extent_yz[3],extent_yz[0],extent_yz[1])
    ax2.imshow(np.flipud(grid_max[x_grid,:,:]),
               extent=extent_zy,cmap=scmap,vmin=lcc_min,vmax=lcc_max,
               rasterized=True)
    ax22.imshow(np.flipud(proj_pdf[x_grid,:,:]),
                extent=extent_zy,cmap=scmap,vmin=lcc_min,vmax=lcc_max,
                rasterized=True)
                      
    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax2.scatter(grid1.z_array[z_max],grid1.y_array[y_max],
                    marker='*', s = 200, linewidths=1,c='g')                      
        ax22.scatter(grid1.z_array[z_max],grid1.y_array[y_max],
                    marker='*', s = 200, linewidths=1,c='g') 
    ax2.axis('tight')
    ax2.set_xlim(Zmin,Zmax)
    ax2.set_ylim(Ymin,Ymax)
    ax2.set_xlabel('Z[km]')
    labels = ax2.get_xticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    labels = ax2.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)

    ax22.axis('tight')
    ax22.set_xlim(Zmin,Zmax)
    ax22.set_ylim(Ymin,Ymax)
    ax22.set_xlabel('Z[km]')
    labels = ax22.get_xticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    labels = ax22.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)

#--plotting vertical projections of the stacked grid
    ax3 = fig.add_axes([0.03, 0.54, 0.4,0.065])
    ax33 = fig.add_axes([0.03, 0.041, 0.4,0.065])
    n_y = np.argwhere(grid1.y_array >= y_eq[0])[0]
    
    ax3.imshow(np.flipud(np.transpose(grid_max[:,y_grid,:])),
               extent=extent_xz,cmap=scmap,vmin=lcc_min,vmax=lcc_max,
               rasterized=True)
    ax33.imshow(np.flipud(np.transpose(proj_pdf[:,y_grid,:])),
                     extent=extent_xz,cmap=scmap,
                      vmin=lcc_min,vmax=lcc_max,
                      rasterized=True)
    
    if grid_max[x_max,y_max,z_max] >= LTrig:
        ax3.scatter(grid1.x_array[x_max],grid1.z_array[z_max],
                    marker='*', s = 200, linewidths=1,c='g')                      
        ax33.scatter(grid1.x_array[x_max],grid1.z_array[z_max],
                    marker='*', s = 200, linewidths=1,c='g')

    ax3.axis('tight')
    ax3.set_xlim(Xmin,Xmax)
    ax3.set_ylim(Zmin,Zmax)
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ####
    #ax3.set_xlabel('X[km]')
    ax3.set_ylabel('Z[km]')

    ax33.axis('tight')
    ax33.set_xlim(Xmin,Xmax)
    ax33.set_ylim(Zmin,Zmax)
    ax33.set_ylim(ax33.get_ylim()[::-1])
    ####
    #ax3.set_xlabel('X[km]')
    ax33.set_ylabel('Z[km]')
    ####

    ax4 = fig.add_axes([0.545,0.44,0.40,0.55])
    for Record,CH_fct,ysta,xsta in zip(st, st_CF,y_sta,x_sta):
        if plot_waveforms:
            ax4.plot(time,(Record.data/max(abs(Record.data)))*10+ysta,'k',alpha=0.4,rasterized=True)
        note1=Record.id
        ax4.plot(time_env,(CH_fct.data/max(abs(CH_fct.data)))*15+ysta,'k',rasterized=True)        
        note2=Record.stats.station
        ax4.text(max(time),ysta,note1,fontsize=10)
        ax1.text(xsta+2,ysta+2, note2,fontsize=12,color='k')
        ax11.text(xsta+2,ysta+2, note2,fontsize=12,color='k')

    ax4.axvspan(t_b,t_e,facecolor='g',alpha=0.2)
    note_t='CF of MBFilter; Fq= '+str(np.round(fq[n22]))+\
            '-'+str(np.round(fq[n1]))+' Hz'
    #fq_str=str(np.round(fq[n1]))+'_'+str(np.round(fq[n22]))
    ax4.text(min(time)+15,65,note_t,fontsize=15)
    
    if len(grid_max[grid_max>=LTrig])>1:
        ax4.axvline(t_b+time_lag/2,linewidth=4, color='r',alpha=0.15)
    
    ax4.set_xlim(0,max(time))
    #ax4.set_xlim(500,1000)
    #ax4.set_ylim(-70,20)
    ax4.set_ylim(-90,80)
    ax4.set_xlabel('Time[sec]')
    ax4.set_ylabel('Y[km]')

    labels = ax4.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)

    file_out_fig = datestr + '_t' +\
                   str('%0.0f' % (t_b)) + 's_' + fq_str + '_fig.png'
    file_out_fig = os.path.join(out_dir, file_out_fig)

    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(file_out_fig)
