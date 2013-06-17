from distutils.core import setup, Extension

module1 = Extension('lib_rec_filter',
                    sources = ['tatka_modules/lib_rec_filter.c'])
module2 = Extension('lib_rec_rms',
                    sources = ['tatka_modules/lib_rec_rms.c'])
module3 = Extension('lib_rec_kurtosis',
                    sources = ['tatka_modules/lib_rec_kurtosis.c'])
module4 = Extension('lib_rec_cc',
                    sources = ['tatka_modules/lib_rec_cc.c'])


setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1,module2,module3,module4])
