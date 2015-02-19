from bp_types import RecursiveMemory

def init_recursive_memory(config):

    n_bands = config.n_freq_bands
    nsamples = int(config.time_lag / config.delta)
    overlap = int(config.t_overlap / config.delta)

    # Create a dictionary of memory objects
    rec_memory = dict()
    for sta_wave in ((sta, wave) for wave in config.wave_type for sta in config.stations):
        # Each entry of the dictionary is a list of memory objects (with n_bands elements)
        rec_memory[sta_wave] = [RecursiveMemory(sta=sta_wave[0], wave=sta_wave[1],
                                                nsamples=nsamples, overlap=overlap)
                                for _ in xrange(n_bands)]

    return rec_memory
