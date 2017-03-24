# -*- coding: utf-8 -*-
"""setup.py: setuptools control."""
from setuptools import setup
from distutils.core import Extension

import inspect
import os
import sys

# Import the version string.
path = os.path.join(os.path.abspath(os.path.dirname(inspect.getfile(
    inspect.currentframe()))), 'backtrackbb')
sys.path.insert(0, path)
from version import get_git_version

with open('README.md', 'rb') as f:
    long_descr = f.read().decode('utf-8')

ext_modules = []
ext_modules.append(Extension(
    'lib_rec_filter',
    sources=['backtrackbb/c_libs/lib_rec_filter.c']))
ext_modules.append(Extension(
    'lib_rec_rms', sources=['backtrackbb/c_libs/lib_rec_rms.c']))
ext_modules.append(Extension(
    'lib_rec_hos', sources=['backtrackbb/c_libs/lib_rec_hos.c']))
ext_modules.append(Extension(
    'lib_rec_cc', sources=['backtrackbb/c_libs/lib_rec_cc.c']))
ext_modules.append(Extension(
    'lib_map_project',
    sources=['backtrackbb/c_libs/map_project/util.c',
             'backtrackbb/c_libs/map_project/map_project.c',
             'backtrackbb/c_libs/map_project/coord_convert.c']))
ext_modules.append(Extension(
    'lib_rosenberger',
    sources=['backtrackbb/c_libs/rosenberger/IA_Kdiag.c',
             'backtrackbb/c_libs/rosenberger/IA_Err.c',
             'backtrackbb/c_libs/rosenberger/IA_Ealloc.c',
             'backtrackbb/c_libs/rosenberger/IA_R2upd.c',
             'backtrackbb/c_libs/rosenberger/rosenberger.c']))

setup(
    name='backtrackbb',
    packages=['backtrackbb', 'backtrackbb.scripts', 'backtrackbb.configobj'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'btbb = backtrackbb.scripts.btbb:main',
            'mbf_plot = backtrackbb.scripts.mbf_plot:main',
            'bt2eventdata = backtrackbb.scripts.bt2eventdata:main',
            'group_triggers = backtrackbb.scripts.group_triggers:main',
            ]
        },
    version=get_git_version(),
    ext_package='backtrackbb.lib',
    ext_modules=ext_modules,
    description='Multi-band array detection and location of seismic sources',
    long_description=long_descr,
    author='Natalia Poiata',
    author_email='poiata@ipgp.fr',
    url='http://backtrackbb.github.io',
    license='CeCILL Free Software License Agreement, Version 2.1',
    platforms='OS Independent',
    classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: CeCILL Free Software License '
                'Agreement, Version 2.1',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics'],
    install_requires=['obspy>=1.0.0', 'scipy>=0.17']
    )
