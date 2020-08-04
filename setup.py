from distutils.core import setup, Extension

pypp = Extension('PyPrinseq',
                    sources = ['src/pyprimer-predictions.c',
                               'src/primer-predictions.c',
                              'src/pyprimer-trimming.c',
                              'src/primer-trimming.c'])

setup (name = 'PyPrinseq',
       version = '1.0',
       description = 'Cleaning fastq sequences. Fast',
       ext_modules = [pypp])
