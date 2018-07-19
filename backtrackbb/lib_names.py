# -*- coding: utf8 -*-
import os
import sysconfig


def get_lib_path(lib):
    suffix = sysconfig.get_config_var('EXT_SUFFIX')
    if not suffix:
        suffix = sysconfig.get_config_var('SO')
    if os.name == 'nt':
        py_version_nodot = sysconfig.get_config_var('py_version_nodot')
        platform = sysconfig.get_platform().replace('-', '_')
        suffix = '.cp{}-{}{}'.format(py_version_nodot, platform, suffix)
    libname = lib + suffix
    return os.path.join(os.path.dirname(__file__), 'lib', libname)
