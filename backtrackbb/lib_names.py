# -*- coding: utf8 -*-
import os
import sysconfig


def get_lib_path(lib):
    suffix = sysconfig.get_config_var('EXT_SUFFIX')
    if not suffix:
        suffix = sysconfig.get_config_var('SO')
    libname = lib + suffix
    return os.path.join(os.path.dirname(__file__), 'lib', libname)
