from distutils.core import setup, Extension

module1 = Extension('lib_rec_filter',
                    sources = ['backtrackbb/c_libs/lib_rec_filter.c'])
module2 = Extension('lib_rec_rms',
                    sources = ['backtrackbb/c_libs/lib_rec_rms.c'])
module3 = Extension('lib_rec_hos',
                    sources = ['backtrackbb/c_libs/lib_rec_hos.c'])
module4 = Extension('lib_rec_cc',
                    sources = ['backtrackbb/c_libs/lib_rec_cc.c'])
module5 = Extension('lib_map_project',
                    sources = ['backtrackbb/c_libs/map_project/util.c',
                               'backtrackbb/c_libs/map_project/map_project.c',
                               'backtrackbb/c_libs/map_project/coord_convert.c'])
module6 = Extension('lib_rosenberger',
                    sources = ['backtrackbb/c_libs/rosenberger/IA_Kdiag.c',
                               'backtrackbb/c_libs/rosenberger/IA_Err.c',
                               'backtrackbb/c_libs/rosenberger/IA_Ealloc.c',
                               'backtrackbb/c_libs/rosenberger/IA_R2upd.c',
                               'backtrackbb/c_libs/rosenberger/rosenberger.c'])


setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1, module2, module3, module4, module5, module6])
