# -*- coding: utf-8 -*-
import setuptools
# Usefull to build f90 files
from numpy.distutils.core import Extension, setup

with open("README.md", "r") as fh:
    long_description = fh.read()

"""
The Fortran extension is build into two steps:
  1. Create a pyf file to avoid trouble with user_defined
fortran type use in the fortran lib. This *.pyf file must be
updated if `wrapper.f90` is modified :
```
f2py -m polsys -h polsys.pyf wrapper.f90 --overwrite-signature 
``
  2. Build the fortran sources and the pyf file,
  using `numpy.distutils.core` and `Extension`.

"""
ext_wrapper = Extension(name='pypolsys.polsys',
                        sources=['pypolsys/src/polsys.pyf',
                                 'pypolsys/801/polsys_plp.f90',
                                 'pypolsys/src/wrapper.f90'],
                        extra_compile_args=['-O3',
                                            '-ffast-math',
                                            '-funroll-loops'],  # "-fopenmp"
                        extra_link_args=['-llapack',
                                         '-lblas'])  # "-lgomp")

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
    ext_modules=[ext_wrapper],
    install_requires=['numpy', 'scipy', 'sympy'],
    classifiers=[
        "Programming Language :: Python :: 3.5",
        "License :: OSI Approved :: GPL 3",
        "Operating System :: OS Independent"],
    # tested with python 3.5 may works with previous py3 version...
    python_requires='>=3.5',
)
