# -*- coding: utf8 -*-
import sys
import os
from configobj import ConfigObj
from validate import Validator
from Config import Config


def _str2bool(arg):
    v = str(arg)
    return v.lower() in ("yes", "true", "t", "1")


def __parse_configspec():
    try:
        configspec_file = os.path.join(os.path.dirname(__file__), 'configspec.conf')
        configspec = ConfigObj(configspec_file, interpolation=False,
                        list_values=False, _inspec=True, file_error=True)
    except IOError, message:
        sys.stderr.write('%s\n' % message)
        sys.exit(1)
    except Exception, message:
        sys.stderr.write('Unable to read "%s": %s\n' % (configspec_file, message))
        sys.exit(1)
    return configspec


def __write_sample_config(configspec, progname):
    c = ConfigObj(configspec=configspec)
    val = Validator()
    c.validate(val)
    c.defaults = []
    c.initial_comment = configspec.initial_comment
    c.comments = configspec.comments
    configfile = progname + '.conf'
    c.write(open(configfile, 'w'))
    print 'Sample config file written to: ' + configfile


def parse_config(config_file):
    configspec = __parse_configspec()
    try:
        config_obj = ConfigObj(config_file, configspec=configspec, file_error=True)
    except IOError, message:
        sys.stderr.write('%s\n' % message)
        sys.exit(1)
    except Exception, message:
        sys.stderr.write('Unable to read "%s": %s\n' % (config_file, message))
        sys.exit(1)

    val = Validator()
    test = config_obj.validate(val)
    if test is True:
        # no problem
        pass
    elif test is False:
        sys.stderr.write('No configuration value present!\n')
        sys.exit(1)
    else:
        for entry in test:
            if test[entry] == False:
                sys.stderr.write('Invalid value for "%s": "%s"\n'
                        % (entry, config_obj[entry]))
        sys.exit(1)

    # Fields needing special treatment:
    if config_obj['save_projGRID'] != 'trigger_only':
        config_obj['save_projGRID'] = _str2bool(config_obj['save_projGRID'])

    stations = config_obj['stations']
    if config_obj['hos_sigma'] is not None:
        hos_sigma = config_obj['hos_sigma']
        # change hos_sigma elements to float
        # and take the power of two
        hos_sigma = [float(x)**2 for x in hos_sigma]
        # make hos_sigma the same length than stations
        if len(stations) > len(hos_sigma):
            hos_sigma += [hos_sigma[-1], ] *\
                (len(stations) - len(hos_sigma))
    else:
        # just create a list of Nones
        hos_sigma = [None, ] * len(stations)
    hos_sigma_dict = {key: value for (key, value) in
        zip(stations, hos_sigma)}
    config_obj['hos_sigma'] = hos_sigma_dict

    # Make wave_type a list
    if config_obj['wave_type'] == 'PS':
        config_obj['wave_type'] = ['P', 'S']
    else:
        config_obj['wave_type'] = [config_obj['wave_type'], ]

    # Create a Config object
    config = Config(config_obj.dict().copy())

    return config
