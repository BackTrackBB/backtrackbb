# -*- coding: utf8 -*-
import os
import importlib.machinery


def get_lib_path(lib):
    suffix = importlib.machinery.EXTENSION_SUFFIXES[0]
    libname = lib + suffix
    return os.path.join(os.path.dirname(__file__), 'lib', libname)
