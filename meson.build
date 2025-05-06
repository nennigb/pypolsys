project(
    'pypolsys',
    ['c', 'fortran'],
    version: run_command('pypolsys/version.py', check: true).stdout().strip(),
    license: 'GNU General Public License v3 (GPLv3)',
    meson_version: '>=0.64.0',
    default_options: ['warning_level=2', 'buildtype=release'],
)


add_languages('fortran', native: true)
# load meson python module
py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()
# https://github.com/numpy/numpy/issues/25000
quadmath_dep = meson.get_compiler('fortran', native: true).find_library('quadmath', required:false)

# Get include from numpy and f2py
# based on https://numpy.org/doc/stable/f2py/buildtools/meson.html
incdir_numpy = run_command(
    py,
    [
        '-c',
        'import os; os.chdir(".."); import numpy; print(numpy.get_include())',
    ],
    check: true,
).stdout().strip()

incdir_f2py = run_command(
    py,
    [
        '-c',
        'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())',
    ],
    check: true,
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

# Need to call f2py to process pyf files
# Based on https://numpy.org/doc/stable/f2py/buildtools/meson.html and
# https://github.com/scipy/scipy/blob/5e986a24f0ab03fbb382ed67ef15895024334481/scipy/io/meson.build#L4
# From https://numpy.org/doc/stable/f2py/buildtools/index.html
# For Fortran 90 programs, f2py generates for an input file blah.f90
#  blahmodule.c
#  blah-f2pywrappers.f, makes subroutines for functions, it rewrites assumed shape arrays as automatic arrays (empty for pypolsys).
#  blah-f2pywrappers2.f90, handles code which is subdivided into modules

# Copy .f2py_f2cmap in buildir (f2py needs it in the working dir)
fs = import('fs')
f2py_f2cmap = fs.copyfile('.f2py_f2cmap')
_polsys_module = custom_target(
    'polsys',
    output: [
        'polsysmodule.c',
        'polsys-f2pywrappers2.f90',
        'polsys-f2pywrappers.f',
    ],
    input: ['pypolsys/src/polsys.pyf'],
    depends: f2py_f2cmap,
    command: [py, '-m', 'numpy.f2py', '@INPUT@', '--lower'],
)
# Pypolsys extension own src files
src_files = [
    _polsys_module,
    'pypolsys/801/polsys_plp.f90',
    'pypolsys/src/wrapper.f90',
]

# Look for lapack only on linux
# Check the operating system the code is runnig (will failed for cross compilation=)
deps = [py_dep, quadmath_dep]
if host_machine.system() == 'linux'
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
      # Need to compile lapack too
      src_files += 'pypolsys/801/lapack_plp.f'
  else
      deps += lapack_dep    
  endif
else
  lapack_dep = disabler()
  # Need to compile lapack too
  src_files += 'pypolsys/801/lapack_plp.f'
endif

# Assumes that f2py was run before in _polsys_module
py.extension_module(
    'polsys',
    src_files,
    incdir_f2py / 'fortranobject.c',
    include_directories: inc_np,
    # No problem if lapack_dep is not found, meson just ignore it
    dependencies: deps,
    install: true,
    native: true,
    subdir: 'pypolsys',  # Folder relative to site-packages to install to
)

# Add pure python files
# Another approach is to split in two meso.build files
python_sources = [
    'pypolsys/__init__.py',
    'pypolsys/test.py',
    'pypolsys/utils.py',
    'pypolsys/version.py',
]

# Install pure Python
# see https://mesonbuild.com/Python-module.html
py.install_sources(
    python_sources,
    pure: false,  # Will be installed next to binaries
    subdir: 'pypolsys',  # Folder relative to site-packages to install to
)

# Add package_data used for tests as source files
install_dir = py.get_install_dir(subdir: 'pypolsys/examples/data')
install_data('pypolsys/examples/data/toy_model.npz', install_dir: install_dir)  # / 'examples' / 'data'

