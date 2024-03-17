project('pypolsys', 'c',
  version : run_command('pypolsys/version.py', check: true).stdout().strip(),
  license: 'BSD-3',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2', 'buildtype=release']
)

add_languages('fortran', native: true)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

# Get include from numpy and f2py
incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

# Need to add f2py pyf files files
# Based on https://numpy.org/doc/stable/f2py/buildtools/meson.html and
# https://github.com/scipy/scipy/blob/5e986a24f0ab03fbb382ed67ef15895024334481/scipy/io/meson.build#L4
# From https://numpy.org/doc/stable/f2py/buildtools/index.html
# For Fortran 90 programs, f2py generates
#        Input file blah.f90
#        Generates:
#            blahmodule.c
#            blah-f2pywrappers.f -> empty for pypolsys
#            blah-f2pywrappers2.f90
# The f90 wrapper is used to handle code which is subdivided into modules. The f wrapper makes subroutines for functions. It rewrites assumed shape arrays as automatic arrays.
_polsys_module = custom_target('polsys',
  output: ['polsysmodule.c', 'polsys-f2pywrappers2.f90', 'polsys-f2pywrappers.f'],
  input: 'pypolsys/src/polsys.pyf',
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '--lower']
)
# Pypolsys own src files
src_files = [_polsys_module,
           'pypolsys/801/polsys_plp.f90',
           'pypolsys/src/wrapper.f90']
# Lapack
lapack_dep = dependency('lapack', required : false)
if not lapack_dep.found()
  # Need to compile lapack too
  src_files += 'pypolsys/801/lapack_plp.f'
endif

# Assumes that f2py was run before in _polsys_module
py.extension_module('polsys', src_files,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  # No problem if lapack_dep is not found, meson just ignore it
  dependencies : [py_dep, lapack_dep],
#  fortran_args: ['-fno-range-check', '-ftree-vectorize', '-funroll-all-loops', '-cpp'],
  install: true,
  link_language: 'fortran',
  native: true
)



