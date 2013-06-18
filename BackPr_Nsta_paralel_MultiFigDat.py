#!/usr/bin/env python
import sys
import os
import numpy as np
import scipy as sp
from obspy.core import read
import itertools
from tatka_modules.mod_filter_picker import make_LinFq, MBfilter_CF
from tatka_modules.NLLGrid import NLLGrid
from tatka_modules.grid_projection import sta_GRD_Proj
from tatka_modules.mod_utils import read_locationTremor,read_locationEQ
from tatka_modules.plot import bp_plot
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
    print "File {0} does not exist".format(inp_file)
    sys.exit(1)

#---Input parameters for BProj run----------------------------------------
ConfigSpec_file = 'configspecBproj.txt'
#Config_file = 'configBproj.txt'
config = Config_Input(ConfigSpec_file,Config_file)

#--input files and station names------------------------------------------
sta=list(set(config.stations))
n_sta = sp.arange(0,len(sta))

t_bb=np.arange(config.start_t,config.end_t,config.t_overlap)
print 'number of time windows=',len(t_bb)

loc_infile = config.catalog_folder+config.data_day+config.tremor_file
location_JMA = config.catalog_folder+config.eq_file
#-------------------------------------------------------------------------       

#--Reading grids of the theoretical travel-times------------------------
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

#--Reading data---------------------------------------------------------
comp = sta[:]
comp[0] = '*'+sta[0]+'*'+config.ch+'*.sac'

st = read(os.path.join(config.data_folder,config.data_day,config.data_hours,comp[0]))
for i in range(1,len(sta)):
    comp[i] = '*'+sta[i]+'*'+config.ch+'*.sac'
    st += read(os.path.join(config.data_folder,config.data_day,config.data_hours,comp[i]))
print 'No of stations in stream = ', len(st)

#-- cut the data to the selected length dt------------------------------
if config.cut_data:
    st.trim(st[0].stats.starttime+0.2*config.delta_t*60.,
			st[0].stats.starttime+1.0*config.delta_t*60.)

    #st.trim(st[0].stats.starttime+0.0*dt,st[0].stats.starttime+1.0*dt)
    #st.trim(st[0].stats.starttime+1.0*dt,st[0].stats.endtime)

#-- decimate to 50Hz sampling-------------------------------------------
for tr in st:
    f_tr = tr.stats.sampling_rate
    if f_tr > config.sr_data:
        dec_ct = int(f_tr/config.sr_data)
        st.decimate(dec_ct, strict_length=False, no_filter=True)

#--remove mean and trend------------------------------------------------
st.detrend(type='constant')
st.detrend(type='linear')

#--Some simple parameters from trace------------------------------------
time = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate
dt=st[0].stats.delta
fs_data = st[0].stats.sampling_rate
dT=dt
npts_d = st[0].stats.npts
n_win_k=config.w_kurt_s /dt
st_CF=st.copy()

#--Calculating frequencies for MBFilter---------------------------------
fq = make_LinFq(config.fq1,config.fq2,npts_d,dT,config.d)
n1=0
n2=len(fq)
n22=len(fq)-1

#----Performing MB filtering and calculating Summary characteristic functions:
for Record,CH_fct in zip(st, st_CF):
    HP2,env_rec,Tn2, Nb2 = MBfilter_CF( Record.data, fq, dT, n_win_k,
                                       CF_type = config.ch_function, var_w = config.win_type )
    CF=env_rec[n1:n2]

    if config.ch_function=='envelope':
        CH_fct.data = np.sqrt((np.sum(CF, axis=0)**2)/len(Tn2[n1:n2]))
    if config.ch_function=='kurtosis':
        CH_fct.data = np.amax(env_rec,axis=0)
        
#----resampling envelopes if desired---------------------------------------------
if config.sr_env < fs_data:
    st_CF.resample(config.sr_env)
    
time_env = np.arange(st_CF[0].stats.npts) / st_CF[0].stats.sampling_rate
dt_env=st_CF[0].stats.delta
npts_e = st_CF[0].stats.npts
fs_env = st_CF[0].stats.sampling_rate
sttime_env = st_CF[0].stats.starttime

print 'frequencies for filtering in (Hz):',fq[n1:n2]

#-------Defining number of station-pairs for calculating LCC---------------------
comb_sta = list(itertools.combinations(sta,2))

#--grid infromation--------------------------------------------------------------
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

x_eq, y_eq, z_eq = read_locationTremor(loc_infile,config.data_hours,
                                       config.lat_orig,config.lon_orig)

##x_jma,y_jma,z_jma = read_locationEQ(location_JMA, config.data_day,config.data_hours,
##                                    config.lat_orig,config.lon_orig)

print 'starting BPmodule'

fq_str=str(np.round(fq[n1]))+'_'+str(np.round(fq[n22]))
#---Loop for moving time window-------------------------------------------------
def run_BackProj(idd):
    t_b=t_bb[idd]
    t_e = t_b+config.time_lag

    stack_grid = np.zeros((nx,ny,nz),float)
    stack_pdf = np.zeros((nx,ny,nz),float)
    
    bb = int(t_b*fs_env)
    ee = int(t_e*fs_env)
    nn = int(config.t_overlap)
    k=0
    #--------------------------------------------------------------------------
    for l in xrange(len(comb_sta)):
        proj_grid = np.zeros((nx,ny,nz),float)
        id1 = sta.index(comb_sta[l][0])
        id2 = sta.index(comb_sta[l][1])
        d = np.sqrt((x_sta[id1]-x_sta[id2])**2+(y_sta[id1]-y_sta[id2])**2)
        if d <= config.maxSTA_distance:
            k+=1
            #print sta[id1],sta[id2],d, l,k
            proj_grid = sta_GRD_Proj(st_CF, GRD_sta,sta, comb_sta, bb, ee, nn,
                                      fs_env,config.time_lag,sttime_env,
                                      config.sm_lcc,nx,ny,nz,l)
            
            stack_grid += proj_grid
            #stack_pdf += 1/np.exp(((1-proj_grid)/proj_grid)**2)
    ## Plotting------------------------------------------------------------------
    bp_plot(grid1, stack_grid/k, comb_sta, x_eq, y_eq,z_eq, config.Trigger,
            t_b, t_e, config.out_dir, config.data_day, config.data_hours, fq_str,
            extent_grd, extent_yz, extent_xz,
            x_sta, y_sta,
            Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
            st, config.scmap, 0.8*config.lcc_max, config.lcc_max,
            sta, st_CF,
            time, time_env, config.time_lag,
            config.plot_waveforms, fq,
            n1, n22)
    
    if config.save_projGRID:
        print 'saving GRID with results'
        out_file = "out_grid/out_"+str(t_b)+".pkl"
        pickle.dump(proj_grid/k, open(out_file, "wb"))

# Running the BackProjection on ncpu--------------------------------------------    
p=Pool(config.ncpu)  #defining number of jobs 
p_outputs = p.map(run_BackProj,xrange(len(t_bb)))
p.close()      #no more tasks
p.join()       #wrap  up current tasks
print 'finished calculations and generating plottings'
