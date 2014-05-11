from distutils.core import setup, Extension

module1 = Extension('lib_rec_filter',
                    sources = ['tatka_modules/c_libs/lib_rec_filter.c'])
module2 = Extension('lib_rec_rms',
                    sources = ['tatka_modules/c_libs/lib_rec_rms.c'])
module3 = Extension('lib_rec_hos',
                    sources = ['tatka_modules/c_libs/lib_rec_hos.c'])
module4 = Extension('lib_rec_cc',
                    sources = ['tatka_modules/c_libs/lib_rec_cc.c'])
module5 = Extension('lib_map_project',
                    sources = ['tatka_modules/c_libs/map_project/util.c',
                               'tatka_modules/c_libs/map_project/map_project.c',
                               'tatka_modules/c_libs/map_project/coord_convert.c'])
module6 = Extension('lib_rosenberger',
                    sources = ['tatka_modules/c_libs/rosenberger/IA_Kdiag.c',
                               'tatka_modules/c_libs/rosenberger/IA_Err.c',
                               'tatka_modules/c_libs/rosenberger/IA_Ealloc.c',
                               'tatka_modules/c_libs/rosenberger/IA_R2upd.c',
                               'tatka_modules/c_libs/rosenberger/rosenberger.c'])


setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1, module2, module3, module4, module5, module6])
