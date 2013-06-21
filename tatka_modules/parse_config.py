# -*- coding: utf8 -*- 
import sys
import os
from configobj import ConfigObj
from validate import Validator
from Config import Config

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
    c.defaults=[]
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

    # Create a Config object
    config = Config(config_obj.dict().copy())

    return config
