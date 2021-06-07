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


# not usefull to map all fortran modules
from .polsys import polsyswrap as polsys
from . import utils
from .utils import solve_univar
from .version import __version__
