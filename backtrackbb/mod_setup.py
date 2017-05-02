# -*- coding: utf8 -*-
"""Setup functions."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import os
from argparse import ArgumentParser
from .configobj import ConfigObj
from .configobj.validate import Validator
from .Config import Config


# Setup ipython shell
if sys.stdout.isatty():
    try:
        import IPython
        ip_version = list(map(int, IPython.__version__.split('.')))
        if ip_version[0] == 0:
            if ip_version[1] >= 11:
                # ipython >= 0.11
                try:
                    from IPython.frontend.terminal.embed \
                        import InteractiveShellEmbed
                    ipshell = InteractiveShellEmbed()
                except ImportError:
                    ipshell = None
            else:
                # ipython < 0.11
                try:
                    from IPython.Shell import IPShellEmbed
                    ipshell = IPShellEmbed()
                except ImportError:
                    ipshell = None
        elif ip_version[0] > 0:
            # ipython >= 1.0.0
            try:
                from IPython.terminal.embed import InteractiveShellEmbed
                ipshell = InteractiveShellEmbed()
            except ImportError:
                ipshell = None
    except ImportError:
        ipshell = None
else:
    ipshell = None


def _parse_args(progname):
    """Parse command line arguments."""
    parser = ArgumentParser()
    parser.add_argument('config_file', metavar='config_file', type=str,
                        help='configuration file')
    if (progname == 'backtrack2eventdata' or
            progname == 'group_triggers'):
        parser.add_argument('trigger_file', metavar='trigger_file', type=str,
                            help='trigger file (output from btbb)')
    if progname == 'backtrack2eventdata':
        parser.add_argument('station_file', metavar='station_file', type=str,
                            help='file with station coordinates (optional)',
                            nargs='?')

    options = parser.parse_known_args()

    return options


def _str2bool(arg):
    v = str(arg)
    return v.lower() in ("yes", "true", "t", "1")


def _parse_configspec():
    try:
        configspec_file = os.path.join(
            os.path.dirname(__file__), 'configspec.conf')
        configspec = ConfigObj(configspec_file, interpolation=False,
                               list_values=False, _inspec=True,
                               file_error=True)
    except IOError as message:
        sys.stderr.write('%s\n' % message)
        sys.exit(1)
    except Exception as message:
        sys.stderr.write('Unable to read "%s": %s\n' %
                         (configspec_file, message))
        sys.exit(1)
    return configspec


def _write_sample_config(configspec, progname):
    c = ConfigObj(configspec=configspec)
    val = Validator()
    c.validate(val)
    c.defaults = []
    c.initial_comment = configspec.initial_comment
    c.comments = configspec.comments
    configfile = progname + '.conf'
    c.write(open(configfile, 'w'))
    print('Sample config file written to: ' + configfile)


def _parse_config(config_file):
    configspec = _parse_configspec()
    try:
        config_obj = ConfigObj(config_file, configspec=configspec,
                               file_error=True)
    except IOError as message:
        sys.stderr.write('%s\n' % message)
        sys.exit(1)
    except Exception as message:
        sys.stderr.write('Unable to read "%s": %s\n' % (config_file, message))
        sys.exit(1)

    # transform strings to one-item lists, when necessary
    for key, val in config_obj.configspec.iteritems():
        if 'list' in val:
            try:
                option = config_obj[key]
            except KeyError:
                continue
            if isinstance(option, str):
                config_obj[key] = [option, ]

    val = Validator()
    test = config_obj.validate(val)
    if isinstance(test, dict):
        for entry in test:
            if not test[entry]:
                sys.stderr.write('Invalid value for "%s": "%s"\n' %
                                 (entry, config_obj[entry]))
        sys.exit(1)
    if not test:
        sys.stderr.write('No configuration value present!\n')
        sys.exit(1)

    # Fields needing special treatment:
    if config_obj['save_projGRID'] != 'trigger_only':
        config_obj['save_projGRID'] = _str2bool(config_obj['save_projGRID'])

    stations = config_obj['stations']
    for hos_sigma_field in ('hos_sigma_P', 'hos_sigma_S'):
        if config_obj[hos_sigma_field] is not None:
            hos_sigma = config_obj[hos_sigma_field]
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
        config_obj[hos_sigma_field] = hos_sigma_dict
    # if there is no value for S, we just make a copy of the dictionary for P
    if all(v is None for v in config_obj['hos_sigma_S'].values()):
        config_obj['hos_sigma_S'] = config_obj['hos_sigma_P'].copy()

    # Make wave_type a list
    if config_obj['wave_type'] == 'PS':
        config_obj['wave_type'] = ['P', 'S']
    else:
        config_obj['wave_type'] = [config_obj['wave_type'], ]

    # Make settings for grid_type
    if config_obj['grid_type'] == 'PS':
        config_obj['grid_type'] = ['P', 'S']
    else:
        config_obj['grid_type'] = config_obj['wave_type']

    # Grid power
    try:
        config_obj['grid_power'] = eval(config_obj['grid_power'],
                                        {'__builtins__': None},
                                        {'nsta': len(config_obj['stations'])})
    except Exception:
        sys.stderr.write('Unable to parse "grid_power". Using 1.\n')
        config_obj['grid_power'] = 1
    if config_obj['grid_power_ellipsoid'] is None:
        config_obj['grid_power_ellipsoid'] = config_obj['grid_power']
    else:
        try:
            config_obj['grid_power_ellipsoid'] =\
                    eval(config_obj['grid_power_ellipsoid'],
                         {'__builtins__': None},
                         {'nsta': len(config_obj['stations'])})
        except Exception:
            sys.stderr.write('Unable to parse "grid_power_ellipsoid". '
                             'Using "grid_power".\n')
            config_obj['grid_power_ellipsoid'] = config_obj['grid_power']
    if config_obj['trigger'] is not None:
        config_obj['trigger'] **= config_obj['grid_power']
    if config_obj['trigger_ellipsoid'] is not None:
        config_obj['trigger_ellipsoid'] **= config_obj['grid_power']

    # Create a Config object
    config = Config(config_obj.dict().copy())

    return config


def configure(progname=None):
    options = _parse_args(progname)
    config = _parse_config(options[0].config_file)
    #add options to config:
    config.options = options[0]
    return config
