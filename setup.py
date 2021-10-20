#!/usr/env python3
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'gfortran decorana_utils.f90 -c -o decorana_utils.o -O3 -fPIC'
print(fortran_mod_comp)
system(fortran_mod_comp)
shared_obj_comp = 'gfortran pydecorana_utils.f90 -c -o pydecorana_utils.o -O3 -fPIC'
print(shared_obj_comp)
system(shared_obj_comp)

ext_modules = [Extension(# module name:
                         'decorana_utils',
                         # source file:
                         ['pydecorana_utils.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['decorana_utils.o', 'pydecorana_utils.o'])]

setup(name = 'decorana_utils',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)