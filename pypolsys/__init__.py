# -*- coding: utf-8 -*-
# This file is part of pypolsys, A python wrapper to `POLSYS_PLP`
# fortran90 package from Layne T. Watson, Steven M. Wise,
# Andrew J. Sommese, August, 1998.

# pypolsys is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# pypolsys is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with pypolsys.  If not, see <https://www.gnu.org/licenses/>.

r"""
pypolsys -- A python wrapper to `POLSYS_PLP` fortran90 homotopy solver
======================================================================
.. include::../README.md
    :start-line:4
    :raw:
"""

# On windows (gfortran+msvc) the fortran dll are put in a `.libs` folder
# by numpy.distutils and need to be in the path before import
# see https://pav.iki.fi/blog/2017-10-08/pywingfortran.html#building-python-wheels-with-fortran-for-windows

import os as _os
import sys as _sys
from platform import system as _system
if _system() == "Windows":
    extra_dll_dir = _os.path.join(_os.path.dirname(__file__), '.libs')
    # For python < 3.8 use path
    if _os.path.isdir(extra_dll_dir):
        _os.environ["PATH"] += _os.pathsep + extra_dll_dir
    # From python >= 3.8 ddl are managed with `os.add_dll_directory`
    if _sys.version_info >= (3, 8):
        if _os.path.isdir(extra_dll_dir):
            _os.add_dll_directory(extra_dll_dir)

# Not usefull to map all fortran modules
from .polsys import polsyswrap as polsys
from . import utils
from .utils import solve_univar
from .version import __version__

# Patch the docstring of the fortran module `polsyswrap` instead of the
# `polsys.xxx.so` file for pdoc3
# Rendering of f2py docstring is not perfect, but better than nothing.
__pdoc__ = {}
__pdoc__['polsys'] = polsys.__doc__
