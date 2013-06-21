#!/usr/bin/env python
import sys
import os
import numpy as np
import itertools
from collections import defaultdict
from tatka_modules.read_traces import read_traces
from tatka_modules.mod_filter_picker import make_LinFq, make_LogFq, MBfilter_CF
from tatka_modules.NLLGrid import NLLGrid
from tatka_modules.mod_utils import read_locationTremor,read_locationEQ
from tatka_modules.grid_projection import sta_GRD_Proj
from tatka_modules.plot import bp_plot, plt_SummaryOut
from tatka_modules.parse_config import parse_config
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
config = parse_config(Config_file)

#---Reading data---------------------------------------------------------
st, stations = read_traces(config)

#------------------------------------------------------------------------
t_bb=np.arange(config.start_t,config.end_t,config.t_overlap)
print 'Number of time windows = ', len(t_bb)

loc_infile = None
location_jma = None
if config.catalog_dir:
    if config.data_day:
        loc_infile = os.path.join(config.catalog_dir, config.data_day+config.tremor_file)
    location_jma = os.path.join(config.catalog_dir, config.eq_file)
#------------------------------------------------------------------------

#--Reading grids of the theoretical travel-times-------------------------
bname =[]
GRD_sta = {}
coord_sta = {}
for station in stations:
    grid_files = '.'.join(('layer', config.wave_type, station, 'time'))
    grid_files = os.path.join(config.grid_dir, grid_files)
    bname.append(grid_files)
    grid = NLLGrid(grid_files)
    coord_sta[station] = (grid.sta_x, grid.sta_y)
    GRD_sta[station] = grid

#--- cut the data to the selected length dt------------------------------
if config.cut_data:
    st.trim(st[0].stats.starttime+0.2*config.delta_t*60.,
                        st[0].stats.starttime+1.0*config.delta_t*60.)

#---- resample data -----------------------------------------------------
if config.sampl_rate_data:
    for tr in st:
        f_tr = tr.stats.sampling_rate
        if f_tr > config.sampl_rate_data:
            dec_ct = int(f_tr/config.sampl_rate_data)
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
n_win_k = config.decay_const /dt

st_CF=st.copy()

#---Calculating frequencies for MBFilter---------------------------------
if config.band_spacing == 'lin':
    fq = make_LinFq(config.f_min, config.f_max, dT, config.n_freq_bands)
elif config.band_spacing == 'log':
    fq = make_LogFq(config.f_min, config.f_max, dT, config.n_freq_bands)
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
if config.sampl_rate_cf:
    if config.sampl_rate_cf < fs_data:
        st_CF.resample(config.sampl_rate_cf)
else:
    config.sampl_rate_cf = fs_data
    
time_env = np.arange(st_CF[0].stats.npts) / st_CF[0].stats.sampling_rate
dt_env=st_CF[0].stats.delta
npts_e = st_CF[0].stats.npts
fs_env = st_CF[0].stats.sampling_rate
sttime_env = st_CF[0].stats.starttime

print 'frequencies for filtering in (Hz):',fq[n1:n2]

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

#----geographical coordinates of the eq's epicenter----------------------
coord_eq = None
if loc_infile:
    coord_eq = read_locationTremor(loc_infile,config.data_hours,
                                       config.lat_orig,config.lon_orig)
coord_jma = None
if location_jma:
    coord_jma = read_locationEQ(location_jma, config.data_day,config.data_hours,
                                    config.lat_orig,config.lon_orig)
#------------------------------------------------------------------------

print 'starting BPmodule'

# Create out_dir, if it doesn't exist
if not os.path.exists(config.out_dir):
    os.mkdir(config.out_dir)

datestr = st[0].stats.starttime.strftime('%y%m%d%H')

fq_str=str(np.round(fq[n1])) + '_' + str(np.round(fq[n22]))
file_out_base = '_'.join((
    datestr,
    str(len(fq)) + 'fq' + fq_str + 'hz',
    str(config.decay_const) + str(config.sampl_rate_cf) + str(config.smooth_lcc) + str(config.t_overlap),
    config.ch_function,
    config.component,
    config.wave_type,
    'trig'+str(config.Trigger)
    ))

file_out_data = file_out_base + '_OUT2.dat'
file_out_data = os.path.join(config.out_dir, file_out_data)

file_out_fig = file_out_base + '_FIG2.png'
file_out_fig = os.path.join(config.out_dir, file_out_fig)


#--------Defining number of station-pairs for calculating LCC------------
comb_sta = list(itertools.combinations(stations, 2))

#------------------------------------------------------------------------
def run_BackProj(idd):
    t_b=t_bb[idd]
    t_e = t_b+config.time_lag

    stack_grid = np.zeros((nx,ny,nz),float)
    #stack_pdf = np.zeros((nx,ny,nz),float)
    
    nn = int(config.t_overlap)
    
    proj_grid = np.zeros((nx,ny,nz),float)

    arrival_times = defaultdict(list)
    k=0
    for sta1, sta2 in comb_sta:
        proj_grid = np.zeros((nx,ny,nz),float)

        x_sta1, y_sta1 = coord_sta[sta1]
        x_sta2, y_sta2 = coord_sta[sta2]
        
        d = np.sqrt((x_sta1-x_sta2)**2+(y_sta1-y_sta2)**2)
        
        if d <= config.maxSTA_distance:
            k+=1

            proj_grid = sta_GRD_Proj(st_CF, GRD_sta, sta1, sta2, t_b, t_e, nn,
                                      fs_env,config.time_lag, sttime_env,
                                      config.smooth_lcc, nx, ny, nz, arrival_times)
            stack_grid += proj_grid
            #stack_pdf += 1/np.exp(((1-proj_grid)/proj_grid)**2)

    ## Plotting------------------------------------------------------------------
    bp_plot(grid1, stack_grid/k, comb_sta, coord_eq, config.Trigger,
            t_b, t_e, config.out_dir, datestr, fq_str,
            extent_grd, extent_yz, extent_xz,
            coord_sta,
            Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
            st, config.scmap, 0.8*config.lcc_max, config.lcc_max,
            stations, st_CF,
            time, time_env, config.time_lag,
            config.plot_waveforms, fq,
            n1, n22)
    
    #Norm_grid= stack_grid/len(comb_sta)
    Norm_grid= stack_grid/k

    Max_NormGrid = np.where(Norm_grid == np.max(Norm_grid))
    x_max = Max_NormGrid[0][0]
    y_max = Max_NormGrid[1][0]
    z_max = Max_NormGrid[2][0]

    if config.save_projGRID:
        print 'saving GRIDS with results'
        out_file = "out_grid/out_"+str(t_b)+".pkl"
        pickle.dump(Norm_grid, open(out_file, "wb"))

    if Norm_grid[x_max,y_max,z_max] >= config.Trigger:
        #for sta in sorted(arrival_times):
        #    print sta, arrival_times[sta]
        
        out = str(grid1.x_array[x_max])+'  '+str(grid1.y_array[y_max])+\
              ' '+str(grid1.z_array[z_max])+'  '+str(t_b)+' '+str(t_e)
        print out
        
        xx_trig = grid1.x_array[x_max]
        yy_trig = grid1.y_array[y_max]
        zz_trig = grid1.z_array[z_max]
        begt_trigWin = t_b
        endt_trigWin = t_e
        centert_trigWin = t_b+config.time_lag/2
        return xx_trig, yy_trig, zz_trig, begt_trigWin, endt_trigWin, centert_trigWin
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

#-plotting output--------------------------------------------------------
plt_SummaryOut(st_CF, st, config.plot_waveforms, config.ch_function, time_env, time,
               coord_sta,
               x_trig, y_trig, z_trig, beg_trigWin, end_trigWin, center_trigWin,t_bb,
               Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, datestr,
               fq[n1],fq[n22],config.time_lag,
               coord_eq, coord_jma, file_out_fig)
