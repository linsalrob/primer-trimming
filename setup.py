from distutils.core import setup, Extension

pypp = Extension('PyPrinseq',
                 sources = [
                     'src/pyprinseq.c',
                     'src/predictprimers.c',
                     'src/trimprimers.c',
                     'src/pyprimer-predictions.c',
                     'src/pyprimer-trimming.c',
                 ])

setup (name = 'PyPrinseq',
       version = '1.0',
       description = 'Cleaning fastq sequences. Fast',
       ext_modules = [pypp])
