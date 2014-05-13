#!/usr/bin/env python
import sys
import os
import numpy as np
import itertools
from collections import defaultdict
from scipy.ndimage.interpolation import zoom
from tatka_modules.bp_types import Trigger
from tatka_modules.read_traces import read_traces
from tatka_modules.mod_filter_picker import make_LinFq, make_LogFq,\
     MBfilter_CF, GaussConv
from tatka_modules.NLLGrid import NLLGrid
from tatka_modules.map_project import get_transform, rect2latlon
from tatka_modules.mod_utils import read_locationTremor, read_locationEQ
from tatka_modules.grid_projection import sta_GRD_Proj
from tatka_modules.plot import bp_plot, plt_SummaryOut
from tatka_modules.parse_config import parse_config
from tatka_modules.mod_bp_TrigOrig_time import TrOrig_time
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
var_twin = config.varWin_stationPair
print 'use of var time window for location:', var_twin
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
bname = []
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
    st.trim(st[0].stats.starttime+config.cut_start,
            st[0].stats.starttime+config.cut_start+config.cut_delta)
else:
    config.cut_start = 0.

#---- resample data -----------------------------------------------------
if config.sampl_rate_data:
    for tr in st:
        f_tr = tr.stats.sampling_rate
        if f_tr > config.sampl_rate_data:
            dec_ct = int(f_tr/config.sampl_rate_data)
            tr.decimate(dec_ct, strict_length=False, no_filter=True)

#---remove mean and trend------------------------------------------------
st.detrend(type='constant')
st.detrend(type='linear')

#---Some simple parameters from trace------------------------------------
time = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate
delta = st[0].stats.delta
fs_data = st[0].stats.sampling_rate
npts_d = st[0].stats.npts
decay_const = config.decay_const
if config.rosenberger_decay_const is not None:
    rosenberger_decay_const = config.rosenberger_decay_const
else:
    rosenberger_decay_const = config.decay_const
sigma_gauss = int(decay_const/delta/4) # so far fixed
                                       # can be added to control file as a separate variable

#---Calculating frequencies for MBFilter---------------------------------
if config.band_spacing == 'lin':
    frequencies = make_LinFq(config.f_min, config.f_max, delta, config.n_freq_bands)
elif config.band_spacing == 'log':
    frequencies = make_LogFq(config.f_min, config.f_max, delta, config.n_freq_bands)
n1 = 0
n2 = len(frequencies)
n22 = len(frequencies) - 1

#----MB filtering and calculating Summary characteristic functions:------
st_CF = st.copy()
for i, station in enumerate(stations):
    st_select = st.select(station=station)
    CF_tr = st_CF.select(station=station)[0]

    print('Creating characteristic function: station No {}/{}'.format(i+1, len(stations)))
    HP2, env_rec, Tn2, Nb2 = MBfilter_CF(st_select, frequencies, var_w=config.win_type,
                                         CF_type=config.ch_function,
                                         CF_decay_win=decay_const,
                                         rosenberger_decay_win=rosenberger_decay_const)

    CF = env_rec[n1:n2]

    if config.ch_function=='envelope':
        CF_tr.data = np.sqrt((np.sum(CF, axis=0)**2)/len(Tn2[n1:n2]))
    if config.ch_function=='kurtosis':
        kurt_argmax = np.amax(env_rec,axis=0)
        CF_tr.data = GaussConv(kurt_argmax, sigma_gauss)

#-----resampling CF if wanted-------------------------------------------
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

print 'frequencies for filtering in (Hz):', frequencies[n1:n2]

#---grid infromation-----------------------------------------------------
grid1 = GRD_sta.values()[0]
nx, ny, nz = np.shape(grid1.array)

#---init map projection--------------------------------------------------
get_transform(grid1.proj_name,
              grid1.orig_lat, grid1.orig_lon,
              grid1.first_std_paral,
              grid1.second_std_paral,
              grid1.map_rot,
              grid1.ellipsoid)

#----geographical coordinates of the eq's epicenter----------------------
coord_eq = None
if loc_infile:
    coord_eq = read_locationTremor(loc_infile, config.data_hours,
                                       config.lat_orig, config.lon_orig)
coord_jma = None
if location_jma:
    coord_jma = read_locationEQ(location_jma, config.data_day, config.data_hours,
                                    config.lat_orig, config.lon_orig)
#------------------------------------------------------------------------

print 'starting BPmodule'

# Create out_dir, if it doesn't exist
if not os.path.exists(config.out_dir):
    os.mkdir(config.out_dir)

datestr = st[0].stats.starttime.strftime('%y%m%d%H')

fq_str = str(np.round(frequencies[n1])) + '_' + str(np.round(frequencies[n22]))
file_out_base = '_'.join((
    datestr,
    str(len(frequencies)) + 'fq' + fq_str + 'hz',
    str(config.decay_const) + str(config.sampl_rate_cf) + str(config.smooth_lcc) + str(config.t_overlap),
    config.ch_function,
    config.channel,
    config.wave_type,
    'trig'+str(config.trigger)
    ))

file_out_data = file_out_base + '_OUT2.dat'
file_out_data = os.path.join(config.out_dir, file_out_data)

file_out_fig = file_out_base + '_FIG2.' + config.plot_format
file_out_fig = os.path.join(config.out_dir, file_out_fig)


#--------Defining number of station-pairs for calculating LCC------------
comb_sta = list(itertools.combinations(stations, 2))
rec_start_time = st[0].stats.starttime
#------------------------------------------------------------------------
max_time_lag = config.time_lag
#------------------------------------------------------------------------
def run_BackProj(idd):
    Mtau = []
    t_b = t_bb[idd]
    t_e = t_b + config.time_lag

    time = np.arange(st[0].stats.npts) / st[0].stats.sampling_rate
    time_env = np.arange(st_CF[0].stats.npts) / st_CF[0].stats.sampling_rate

    stack_grid = np.zeros((nx,ny,nz),float)
    proj_grid = np.zeros((nx,ny,nz),float)
    arrival_times = defaultdict(list)
    trig_time = defaultdict(list)
    bp_trig_time = defaultdict(list)

    nn = int(config.t_overlap)

    k=0
    for sta1, sta2 in comb_sta:
        proj_grid = np.zeros((nx,ny,nz),float)

        x_sta1, y_sta1 = coord_sta[sta1]
        x_sta2, y_sta2 = coord_sta[sta2]

        distance = np.sqrt((x_sta1-x_sta2)**2+(y_sta1-y_sta2)**2)

        if distance <= config.maxSTA_distance:
            k+=1
            if var_twin:
                tau_max = GRD_sta[sta1].get_value(GRD_sta[sta2].sta_x,
                                                  GRD_sta[sta2].sta_y,
                                                  GRD_sta[sta2].sta_z)
                Mtau.append(np.round(tau_max,1))
                t_e = t_b + np.round(tau_max,1)
            else:
                Mtau = 'none'
                tau_max = None

            proj_grid = sta_GRD_Proj(st_CF, GRD_sta, sta1, sta2, t_b, t_e, nn,
                                      fs_env, sttime_env, config,
                                      nx, ny, nz, arrival_times, tau_max)
            stack_grid += proj_grid

    Norm_grid = stack_grid/k

    Max_NormGrid = np.where(Norm_grid == np.max(Norm_grid))
    i_max = Max_NormGrid[0][0]
    j_max = Max_NormGrid[1][0]
    k_max = Max_NormGrid[2][0]

    if config.cut_data:
        start_tw = config.cut_start+t_b
    else:
        start_tw = t_b
        config.cut_start = 0.

    if config.save_projGRID:
        print 'saving GRIDS with results'
        out_file = "out_grid/out_"+str(t_b)+".pkl"
        pickle.dump(Norm_grid, open(out_file, "wb"))

    trigger = None
    if Norm_grid[i_max, j_max, k_max] >= config.trigger:
        if config.max_subdivide is not None:
            #We "zoom" the grid in a region close to the maximum
            zoom_factor = float(config.max_subdivide)
            #TODO: slicing will not work if the maximum
            #is close to the grid edge!
            slice_factor = 4 #this is hardcoded, for now
            sf = int(slice_factor)
            slice_grid = Norm_grid[i_max-sf:i_max+sf,
                                   j_max-sf:j_max+sf,
                                   k_max-sf:k_max+sf]
            zoom_slice_grid = zoom(slice_grid, zoom_factor)
            zoom_i_max, zoom_j_max, zoom_k_max =\
                    [a[0]/zoom_factor
                     for a in np.where(zoom_slice_grid == np.max(zoom_slice_grid))]
            i_max = zoom_i_max + i_max-sf
            j_max = zoom_j_max + j_max-sf
            k_max = zoom_k_max + k_max-sf
        xx_trig, yy_trig, zz_trig = grid1.get_xyz(i_max, j_max, k_max)

        trigger = Trigger()
        trigger.x, trigger.y, trigger.z = grid1.get_xyz(i_max, j_max, k_max)
        trigger.i, trigger.j, trigger.k = i_max, j_max, k_max
        trigger.max_grid = np.round(np.max(Norm_grid),4)
        trigger.beg_win = start_tw
        trigger.end_win = start_tw + max_time_lag
        trigger.center_win = start_tw + max_time_lag/2.
        trigger.lat, trigger.lon =\
                rect2latlon(trigger.x, trigger.y)

    ##------------------Origin time calculation------------------------------------------------------
        bp_origin_time, bp_trig_time =\
                TrOrig_time(config, stations, GRD_sta, xx_trig, yy_trig, zz_trig,
                            rec_start_time, arrival_times, trig_time)

        trigger.origin_time = bp_origin_time
    ##-----------------------------------------------------------------------------------------------
        print trigger

    ## Plotting------------------------------------------------------------------
    bp_plot(config, grid1, Norm_grid, comb_sta,
            coord_eq, t_b, t_e, datestr, fq_str,
            coord_sta, st, stations, st_CF,
            time, time_env,
            frequencies, n1, n22, trigger, arrival_times, bp_trig_time, Mtau)

    return trigger

#------end loop for BackProj---------------------------------------------

#---running program------------------------------------------------------
if config.ncpu > 1:
    # parallel execution
    print 'Running on %d threads' % config.ncpu
    p = Pool(config.ncpu)  #defining number of jobs
    p_outputs = p.map(run_BackProj,xrange(len(t_bb)))
    p.close()      #no more tasks
    p.join()       #wrap  up current tasks
else:
    # serial execution (useful for debugging)
    print 'Running on 1 thread'
    p_outputs = []
    for idd in xrange(len(t_bb)):
        p_outputs.append(run_BackProj(idd))

triggers = filter(None, p_outputs)

#----------Outputs-------------------------------------------------------
#writing output
with open(file_out_data,'w') as f:
    for trigger in triggers:
        f.write(str(trigger) + '\n')

#-plotting output--------------------------------------------------------
plt_SummaryOut(config, grid1, st_CF, st, time_env, time, coord_sta,
               triggers, t_bb, datestr, frequencies[n1], frequencies[n22],
               coord_eq, coord_jma, file_out_fig)
