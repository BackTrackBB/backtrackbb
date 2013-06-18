#!/usr/bin/env python
import sys
import os
import numpy as np
import scipy as sp
from obspy.core import read
import itertools
from tatka_modules.mod_filter_picker import make_LinFq, MBfilter_CF
from tatka_modules.NLLGrid import NLLGrid
from tatka_modules.mod_utils import read_locationTremor,read_locationEQ
from tatka_modules.grid_projection import sta_GRD_Proj
from tatka_modules.plot import plt_SummaryOut
from tatka_modules.Config import Config_Input
import cPickle as pickle
from multiprocessing import Pool

#-------------------------------------------------------------------------
if len(sys.argv) != 2:
    print "this_code  <input_config_file>"
    sys.exit(1)
else:
    Config_file = sys.argv[1]
if not os.path.isfile(Config_file):
    print "File {0} does not exist".format(Config_file)
    sys.exit(1)

#---Input parameters for BProj run----------------------------------------
ConfigSpec_file = 'configspecBproj.txt'
config = Config_Input(ConfigSpec_file,Config_file)
#------------------------------------------------------------------------
sta=list(set(config.stations))
n_sta = sp.arange(0,len(sta))

t_bb=np.arange(config.start_t,config.end_t,config.t_overlap)
print 'number of time windows=',len(t_bb)

loc_infile = config.catalog_folder+config.data_day+config.tremor_file
location_JMA = config.catalog_folder+config.eq_file
#------------------------------------------------------------------------

#--Reading grids of the theoretical travel-times-------------------------
x_sta=[]
y_sta=[]
bname =[]
GRD_sta=[]
for station in sta:
    grid_files = config.grid_dir+config.wave_type+station+'.time'
    bname.append(grid_files)
    x_sta.append(NLLGrid(grid_files).sta_x)
    y_sta.append(NLLGrid(grid_files).sta_y)
    GRD_sta.append(NLLGrid(grid_files).array)

x_sta=tuple(x_sta)
y_sta=tuple(y_sta)
GRD_sta=tuple(GRD_sta)

#---Reading data---------------------------------------------------------
comp = sta[:]
comp[0] = '*'+sta[0]+'*'+config.ch+'*.'+config.data_type

st = read(os.path.join(config.data_folder,config.data_day,config.data_hours,comp[0]))

for i in xrange(1,len(sta)):
    comp[i] = '*' + sta[i] + '*' + config.ch + '*.'+config.data_type
    st += read(os.path.join(config.data_folder,config.data_day,config.data_hours,comp[i]))
print 'No of stations in stream = ', len(st)

#--- cut the data to the selected length dt------------------------------
if config.cut_data:
    st.trim(st[0].stats.starttime+0.2*config.delta_t*60.,
                        st[0].stats.starttime+1.0*config.delta_t*60.)

#----decimate to 50Hz sampling-------------------------------------------
for tr in st:
    f_tr = tr.stats.sampling_rate
    if f_tr > config.sr_data:
        dec_ct = int(f_tr/config.sr_data)
        st.decimate(dec_ct, strict_length=False, no_filter=True)

#---remove mean and trend------------------------------------------------
st.detrend(type='constant')
st.detrend(type='linear')
#---Some simple parameters from trace------------------------------------
time = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate
dt=st[0].stats.delta
fs_data = st[0].stats.sampling_rate
dT=dt
npts_d = st[0].stats.npts
n_win_k=config.w_kurt_s /dt

st_CF=st.copy()

#---Calculating frequencies for MBFilter---------------------------------
fq = make_LinFq(config.fq1,config.fq2,npts_d,dT,config.d)
n1=0
n2=len(fq)
n22=len(fq)-1

#----MB filtering and calculating Summary characteristic functions:------
for Record,CH_fct in zip(st, st_CF):
    HP2,env_rec,Tn2, Nb2 = MBfilter_CF( Record.data, fq, dT, n_win_k,
                                       CF_type = config.ch_function, var_w = config.win_type )
    CF=env_rec[n1:n2]

    if config.ch_function=='envelope':
        CH_fct.data = np.sqrt((np.sum(CF, axis=0)**2)/len(Tn2[n1:n2]))
    if config.ch_function=='kurtosis':
        CH_fct.data = np.amax(env_rec,axis=0)
            
#-----resampling envelopes if wanted------------------------------------
if config.sr_env < fs_data:
    st_CF.resample(config.sr_env)
    
time_env = np.arange(st_CF[0].stats.npts) / st_CF[0].stats.sampling_rate
dt_env=st_CF[0].stats.delta
npts_e = st_CF[0].stats.npts
fs_env = st_CF[0].stats.sampling_rate
sttime_env = st_CF[0].stats.starttime

print 'frequencies for filtering in (Hz):',fq[n1:n2]

#--------Defining number of station-pairs for calculating LCC------------
comb_sta = list(itertools.combinations(sta,2))

#---grid infromation-----------------------------------------------------
grid1 = NLLGrid(bname[1])
grid1.init_xyz_arrays()
extent_grd = grid1.get_xy_extent()
extent_yz = grid1.get_yz_extent()
extent_xz = (extent_grd[0],extent_grd[1],extent_yz[2],extent_yz[3])
Xmax = max(grid1.x_array)
Xmin = min(grid1.x_array)
Ymax = max(grid1.y_array)
Ymin = min(grid1.y_array)
Zmax = max(grid1.z_array)
Zmin = min(grid1.z_array)
nx,ny,nz = np.shape(grid1.array)

print 'starting BPmodule'
fq_str=str(np.round(fq[n1]))+'_'+str(np.round(fq[n22]))

file_out_data = config.out_dir+config.data_day+config.data_hours+'_'+str(len(fq))+'fq'+\
                fq_str+'hz_'+str(config.w_kurt_s)+str(config.sr_env)+str(config.sm_lcc)+\
                str(config.t_overlap)+'_'+config.ch_function+'_'+config.ch+config.wave_type+\
                'trig'+str(config.Trigger)+'_OUT2.dat'

file_out_fig  = config.out_dir+config.data_day+config.data_hours+'_'+str(len(fq))+'fq'+\
                fq_str+'hz_'+str(config.w_kurt_s)+str(config.sr_env)+str(config.sm_lcc)+\
                str(config.t_overlap)+'_'+config.ch_function+'_'+config.ch+config.wave_type+\
                'trig'+str(config.Trigger)+'_FIG2.png'


#------------------------------------------------------------------------
def run_BackProj(idd):
    t_b=t_bb[idd]
    t_e = t_b+config.time_lag
    
    bb = int(t_b*fs_env)
    ee = int(t_e*fs_env)
    nn = int(config.t_overlap)
    
    proj_grid = np.zeros((nx,ny,nz),float)
    k=0
    for l in xrange(len(comb_sta)):
        id1 = sta.index(comb_sta[l][0])
        id2 = sta.index(comb_sta[l][1])
        
        d = np.sqrt((x_sta[id1]-x_sta[id2])**2+(y_sta[id1]-y_sta[id2])**2)
        
        if d <= config.maxSTA_distance:
            k+=1

            proj_grid += sta_GRD_Proj(st_CF, GRD_sta,sta, comb_sta, bb, ee, nn,
                                      fs_env,config.time_lag,sttime_env,
                                      config.sm_lcc,nx,ny,nz,l)
        
    #Norm_grid= proj_grid/len(comb_sta)
    Norm_grid= proj_grid/k

    Max_NormGrid =np.where(Norm_grid == np.max(Norm_grid))
    x_max = Max_NormGrid[0][0]
    y_max = Max_NormGrid[1][0]
    z_max = Max_NormGrid[2][0]

    if config.save_projGRID:
        print 'saving GRIDS with results'
        out_file = "out_grid/out_"+str(t_b)+".pkl"
        pickle.dump(Norm_grid, open(out_file, "wb"))

    if Norm_grid[x_max,y_max,z_max] >= config.Trigger:
        
        out = str(grid1.x_array[x_max])+'  '+str(grid1.y_array[y_max])+\
              ' '+str(grid1.z_array[z_max])+'  '+str(t_b)+' '+str(t_e)
        print out
        
        xx_trig =grid1.x_array[x_max]
        yy_trig = grid1.y_array[y_max]
        zz_trig = grid1.z_array[z_max]
        begt_trigWin = t_b
        endt_trigWin = t_e
        centert_trigWin = t_b+config.time_lag/2
        return xx_trig, yy_trig,zz_trig,begt_trigWin,endt_trigWin,centert_trigWin
#------end loop for BackProj---------------------------------------------

x_trig=[]
y_trig=[]
z_trig=[]
beg_trigWin=[]
end_trigWin=[]
center_trigWin=[]

#---running program------------------------------------------------------
p=Pool(config.ncpu)  #defining number of jobs 
p_outputs = p.map(run_BackProj,xrange(len(t_bb)))
p.close()      #no more tasks
p.join()       #wrap  up current tasks

f_output = filter(None,p_outputs)

for out in f_output:
    x_trig.append(out[0])
    y_trig.append(out[1])
    z_trig.append(out[2])
    beg_trigWin.append(out[3])
    end_trigWin.append(out[4])
    center_trigWin.append(out[5])
    
#----------Outputs-------------------------------------------------------
#writting output
f = open(file_out_data,'w')
[f.write(str(x_trig[k])+' '+str(y_trig[k])+' '+str(z_trig[k])+' '+str(beg_trigWin[k])+' '+\
         str(end_trigWin[k])+'\n')for k in xrange(len(beg_trigWin))]
f.close()

#----geographical coordinates of the eq's epicenter----------------------
x_eq, y_eq, z_eq = read_locationTremor(loc_infile,config.data_hours,
                                       config.lat_orig,config.lon_orig)

x_jma,y_jma,z_jma = read_locationEQ(location_JMA, config.data_day,config.data_hours,
                                    config.lat_orig,config.lon_orig)
#------------------------------------------------------------------------
#-plotting output--------------------------------------------------------
plt_SummaryOut(st_CF, st, config.plot_waveforms, config.ch_function, time_env, time,
                           sta, x_sta, y_sta,
               x_trig, y_trig, z_trig, beg_trigWin, end_trigWin, center_trigWin,t_bb,
               Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, config.data_day, config.data_hours,
               fq[n1],fq[n22],config.time_lag,
               x_eq, y_eq, z_eq, x_jma, y_jma, z_jma, file_out_fig)
