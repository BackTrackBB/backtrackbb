import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pylab
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

def bp_plot(fig, grid1, grid_max, proj_pdf, x_eq, y_eq,z_eq, LTrig,
        t_b, t_e, out_dir, data_day, data_hours, fq_str,
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
    cb1=plt.colorbar(hnd, cax=cbx1,orientation='horizontal')
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
    cb11=plt.colorbar(hnd2, cax=cbx11,orientation='horizontal')
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

##    file_out_fig  = out_dir+data_day+data_hours+'_t'+\
##                str('%0.0f' % (t_b)) +'s_'+ fq_str+'_fig.png'
##    fig.savefig(file_out_fig, dpi=100)

    file_out_fig  = out_dir+data_day+data_hours+'_t'+\
                str('%0.0f' % (t_b)) +'s_'+ fq_str+'_fig.png'
    fig.savefig(file_out_fig)
    
    plt.close('all')

    
def plt_SummaryOut(fig, st_CF, st, plot_waveforms, ch_function,time_env, time, sta, x_sta, y_sta,
                   x_trig, y_trig, z_trig, beg_trigWin, end_trigWin, center_trigWin,t_bb,
                   Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, data_day, data_hours,fq_1,fq_2,time_lag,
                   x_eq, y_eq, z_eq, x_jma, y_jma, z_jma, file_out_fig):
    
    ax4 = fig.add_axes([0.05,0.472,0.85,0.52])

    for Record,CH_fct,ysta in zip(st, st_CF,y_sta):
        ax4.plot(time_env,(CH_fct.data/max(abs(CH_fct.data)))*10+ysta,
                 'k',rasterized=True)

        if plot_waveforms:
            ax4.plot(time,(Record.data/max(abs(Record.data)))*10+ysta,'k',
                     alpha=0.4,rasterized=True)

        note1=Record.id
        note2=Record.stats.station
        ax4.text(max(time),ysta,note1,fontsize=10)

    note=ch_function+' of MBFilter; Fq. range: '+str(np.round(fq_1))+\
            '-'+str(np.round(fq_2))+' Hz'
    ax4.text(min(time)+15,70,note,fontsize=15)
    ##[ax4.axvline(center_trigWin[m],linewidth=2, color='r',alpha=0.2)\
    ##     for m in xrange(len(beg_trigWin))]

    [ax4.axvspan(beg_trigWin[m],end_trigWin[m],facecolor='r',alpha=0.1)\
         for m in xrange(len(beg_trigWin))]

    ax4.axvline(t_bb[0],linewidth=1, color='b',alpha=0.9)
    ax4.axvline(t_bb[-1]+time_lag,linewidth=1, color='b',alpha=0.9)

    ax4.set_xlim(0,max(time))
    #ax4.set_ylim(-70,20)
    ax4.set_ylim(-90,80)
    ax4.set_xlabel('Time[sec]')
    ax4.set_ylabel('Y[km]')
    labels = ax4.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)


    ax5 = fig.add_axes([0.22,0.19,0.4, 0.25])
    ax5.scatter(x_trig,y_trig,marker='*', s = 80, linewidths=0.5,c='g',alpha=0.7)    
    ax5.scatter(x_sta,y_sta,marker='^', s = 250, linewidths=1,c='k',alpha=0.79)
    ax5.scatter(x_eq,y_eq,marker='*', s = 300, linewidths=1,c='r')
    ax5.scatter(x_jma,y_jma,marker='o', s = 90, linewidths=1,c='m')
    ax5.axis('tight')
    ax5.set_xlim(Xmin,Xmax)
    ax5.set_ylim(Ymin,Ymax)
    ax5.set_xlabel('X[km]')
    ax5.set_ylabel('Y[km]')
    labels = ax5.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)
    note='Day: '+data_day+',  Hour: '+data_hours
    ax5.text(Xmin,Ymax,note,fontsize=15)
    #ax5.set_aspect('equal')

    ax6 = fig.add_axes([0.605,0.19,0.19, 0.25])
    ax6.scatter(z_trig,y_trig,marker='*', s = 80, linewidths=0.5,c='g',alpha=0.7)
    ax6.axis('tight')
    ax6.scatter(z_eq,y_eq,marker='*', s = 300, linewidths=1,c='r')
    ax6.scatter(z_jma,y_jma,marker='o', s = 90, linewidths=1,c='m')
    ax6.set_xlim(Zmin,Zmax)
    ax6.set_ylim(Ymin,Ymax)
    ax6.set_aspect('equal')

    ax6.set_xlabel('Z[km]')
    labels = ax6.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)


    ax7 = fig.add_axes([0.22,0.027,0.4, 0.15])
    ax7.scatter(x_trig,z_trig,marker='*', s = 80, linewidths=0.5,c='g',alpha=0.7)
    ax7.axis('tight')
    ax7.scatter(x_eq,z_eq,marker='*', s = 300, linewidths=1,c='r')
    ax7.scatter(x_jma,z_jma,marker='o', s = 90, linewidths=1,c='m')
    ax7.set_xlim(Xmin,Xmax)
    ax7.set_ylim(Zmin,Zmax)
    ax7.set_ylim(ax7.get_ylim()[::-1])
    ax7.set_aspect('equal')
    ax7.set_xlabel('X[km]')
    ax7.set_ylabel('Z[km]')
    labels = ax7.get_yticklabels()
    pylab.setp(labels, rotation=90, fontsize=12)

    fig.savefig(file_out_fig)
    plt.close('All')
