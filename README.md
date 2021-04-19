This program contains a **minimal** python wrapper to `POLSYS_PLP` fortran90 package from 
Layne T. Watson, Steven M. Wise, Andrew
J. Sommese, August, 1998.


`POLSYS_PLP` is a solver for N complex coefficient polynomial systems of equations in N unknowns by a probability-one, globally convergent homotopy method.

For a quick survey of [homotopy continuation](https://en.wikipedia.org/wiki/Numerical_algebraic_geometry#Homotopy_continuation) or the very clear pratical presentation [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/guides) website.

There several other homotopy softwares (more recent, also with python interface, ...) but such tool will be hard to automatically install and deploy for some applications. 
The advantage of `POLSYS_PLP` is a self content single f90 file. For more advances usages :
  - [PHCpack](https://github.com/janverschelde/PHCpack)
  - [Bertini](https://github.com/bertiniteam/b2)
  - [Hom4PSpy](http://www.hom4ps3.org) and the wrapper on [pypi](https://pypi.org/project/Hom4PSpy/)
  - [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org)
  
# POLSYS_PLP
The whole algorithm is described in
> Wise, Steven M., Andrew J. Sommese, and Layne T. Watson. "Algorithm 801: POLSYS_PLP: A partitioned linear product homotopy code for solving polynomial systems of equations." ACM Transactions on Mathematical Software (TOMS) 26.1 (2000): 176-200.

  - The sources are availlable on [netlib](http://www.netlib.org/opt/polsysplp.tgz).
  - This code is based HOMPACK, the paper is availlable [here](https://vtechworks.lib.vt.edu/bitstream/handle/10919/19612/TR-90-36.pdf)
  - The theoritical part of `POLSYS_PLP` is availlable [here](https://vtechworks.lib.vt.edu/bitstream/handle/10919/36933/thesis.pdf)
  
To make easier the build of this module, a copy of `POLSYS_PLP` *original* sources is included in the `801` folder.
  
# Building
```
python3 -m numpy.f2py -c polsys_plp.f90 -m polsys
```
doesn't work out of the box, but work if I used 
```
gfortran -c -fPIC polsys_plp.f90 -llapack -lblas
f2py -c --fcompiler=gfortran polsys_plp.o -llapack -lblas -m polsys calltest.f90 
```
I found some problem related to some `gfortran` versoin (see [here](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=84276)) and to compile on ubuntu 16.04, I remove `intent(in, out)` statement in the main file.

It seems to work out of the box with `gfortan-8`.


## License
This file is part of pypolsys, a simple python wrapper to fortran package polsys_plp.
pypolsys is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
pypolsys is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with amc2moodle.  If not, see <https://www.gnu.org/licenses/>.
