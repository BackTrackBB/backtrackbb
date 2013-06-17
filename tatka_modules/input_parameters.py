

class InParameters(object):
    '''
    Describes the class of input parameters 
    '''
    def __init__(self):
        self.inpara = {}

        self.inpara['sta'] = None


#####--input files and station names--------------------------------------------
#####sta = ['MISH']
####sta = ['MISH','TSMH','OOZH']
######sta=['MISH','IKTH','OOZH','TBEH','KWBH','YNDH','HIYH','TSYH','TSMH','UWAH',\
######     'NAKH']
####
####n_sta = sp.arange(0,len(sta))
####
####data_folder = 'Tremor_data/data'
####ch = 'E'
####data_day = "120527"
####data_hours = "21"
####
####grid_dir = 'nll_Shikoku/time/layer.S.'
####
####out_dir = ''
####
#####--checking for repeating stations--------------------------------------
####d=[x for x, y in collections.Counter(sta).items() if y > 1]
####if len(d)>0:
####    error_mes='there are duplicate stations in station list: '+str(d)
####    sys.exit(error_mes)
######else:
######    print 'station list check ok'
#####-----------------------------------------------------------------------
####    
####time_lag = 26.                              #defines max_lag in seconds
####w_kurt_s=3.0                                #value of attenuation constant for rec. kurt. in seconds
####ch_function = 'envelope'
######ch_function = 'kurtosis'
####
#####--Parameters for the sliding window------------------------------------
######start_t = 3.
######end_t = 1160.
######end_t = 9.
####
####start_t = 870.
####end_t = 876.
####
####t_overlap = 3.                                #shift of each following time window in calc
####
####t_bb=np.arange(start_t,end_t,t_overlap)
####print 'number of time windows=',len(t_bb)
####
####sr_data = 50.
####
######sr_env = 1
######sm_lcc = 1.2                                #smoothing parameter for LCC    
####sr_env = 50
####sm_lcc = 40                                   #smoothing parameter for LCC
####
#####--Parameters for the MBF analysis--------------------------------------
####Tf=1.                                       # min filter frequency for MB filter
####f_max = 50.
####fq1 = 2.0
####fq2 = 10.
######d = 4000
####d = 1200
####
####
####cut_data = True
####dt = 25.*60.                                #time length of the data to be used in calculations
####
####LTrig = 0.85
####lcc_max =0.95
####lcc_min =0.8*lcc_max
####
####scmap = 'jet'
####
#####--geogr. coordinates for the origin of the coordinate system-----------
####lat_orig = 33.618
####lon_orig = 132.988
#####--geographical coordinates of the eq's epicenter
####lat_eq = 33.3786
####lon_eq = 132.4777
####zh = 30
####
#####--parameters for plotting-----------------------------------------------
####plot_waveforms = True
#####------------------------------------------------------------------------

#------------------------------------------------------------------------
