# -*- coding: utf-8 -*-
"""
# pypolsys setup.

Define the fortran extension and its dependancies.
pypolsys can be either:
    - linked to lapack/blas (default on linux)
    - build with `POLSYS_PLP` lapack sources (default on windows)
The behavior can be specified with the install/developp flag `--buildLapack`

"""
import setuptools
import platform
import sys
# Usefull to build f90 files
from numpy.distutils.core import Extension, setup
# To setup user options
from setuptools.command.install import install
from setuptools.command.develop import develop

class UserOptionCommand(object):
    """ Mixin class to add user option support.

    based on https://stackoverflow.com/questions/18725137/how-to-obtain-arguments-passed-to-setup-py-from-pip-with-install-option
    """
    # Other options may be added here
    user_options = [
        ('buildLapack', None, 'Provide the required lapack/blas sources and build it.')]

    def initialize_options(self):
        """Initialize default options."""
        super().initialize_options()
        self.buildLapack = None

    def finalize_options(self):
        """Validate options."""
        # Nothing to do here
        super().finalize_options()

    def run(self):
        """Use options."""
        print('> Building with lapackBuild =', buildLapack)
        super().run()


# Update setuptool classes
class InstallCommand(UserOptionCommand, install):
    user_options = getattr(install, 'user_options', []) + UserOptionCommand.user_options


class DevelopCommand(UserOptionCommand, develop):
    user_options = getattr(develop, 'user_options', []) + UserOptionCommand.user_options


def ext_wrapper(buildLapack):
    """ Configure the fortran POLSYS_PLP extension.
    """
    # Base files
    sources = ['pypolsys/src/polsys.pyf',
               'pypolsys/801/polsys_plp.f90',
               'pypolsys/src/wrapper.f90']

    # Configure the POLSYS_PLP Extension building
    if buildLapack == 1:
        # Add lapack_plp to sources
        sources.append('pypolsys/801/lapack_plp.f')
        extra_link_args = []
    else:
        # Link to system lapack/blas
        extra_link_args = ['-llapack', '-lblas']
    return Extension(name='pypolsys.polsys',
                     sources=sources,
                     extra_compile_args=['-O3',
                                         '-ffast-math',
                                         '-funroll-loops'],  # "-fopenmp"
                     extra_link_args=extra_link_args)         # "-lgomp")

# Load README
with open("README.md", "r") as fh:
    long_description = fh.read()

# Check platerform and options
buildLapack = None
if '--buildLapack' in sys.argv:
    buildLapack = 1
# On windows, by default, we used the lapack version shipped with POLSYS_PLP
this_os = platform.system().lower()
if this_os == 'windows':
    buildLapack = 1
print('> Building plateform =', this_os)
print('> Building with lapackBuild =', buildLapack)

# Fill setuptools
setup(
    name="pypolsys",
    version="0.1",
    author="B. Nennig",
    author_email="benoit.nennig@supmeca.fr",
    description=r"A python wrapper to `POLSYS_PLP` fortran90 package from Layne T. Watson,\
    Steven M. Wise, Andrew J. Sommese, August, 1998.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nennigb/pypolsys",
    # Use find_packages() to automatically discover all packages and subpackages
    packages=setuptools.find_packages(),
    include_package_data=True,
    # build f90 module
    ext_modules=[ext_wrapper(buildLapack)],
    install_requires=['numpy', 'scipy', 'sympy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL 3",
        "Operating System :: OS Independent"],
    # tested with python 3.5 may works with previous py3 version...
    python_requires='>=3.5',
    # pass modified install and develop class to support user options
    cmdclass={
        'install': InstallCommand,
        'develop': DevelopCommand,
        }
)
