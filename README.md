PYPOLSYS
========
![CI-Ubuntu](https://github.com/nennigb/pypolsys/workflows/CI-Ubuntu/badge.svg) ![CI-Windows](https://github.com/nennigb/pypolsys/workflows/CI-Windows/badge.svg)

This package provides a python wrapper to `POLSYS_PLP` fortran90 package from Layne T. Watson, Steven M. Wise, Andrew J. Sommese, August, 1998.

`POLSYS_PLP` is a solver for N complex coefficients polynomial systems of equations in N unknowns by a probability-one, globally convergent homotopy method.

For a quick survey on [homotopy continuation](https://en.wikipedia.org/wiki/Numerical_algebraic_geometry#Homotopy_continuation) see the very clear practical presentation available on [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/guides) website or this nice [tutorial](http://homepages.math.uic.edu/~jan/tutorial.pdf).

There are several other homotopy softwares (more recent, also with python interface, ...) allowing more advanced applications. See for instance:
  - [PHCpack](https://github.com/janverschelde/PHCpack)
  - [Bertini](https://github.com/bertiniteam/b2)
  - [Hom4PSpy](http://www.hom4ps3.org) and the wrapper on [pypi](https://pypi.org/project/Hom4PSpy/)
  - [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org)

The advantage of `POLSYS_PLP` is to be an open-source self-content single f90 file (it requires few `blas` and `lapack` functions that are shipped with `POLSYS_PLP` sources or that may be linked to optimized libraries). The idea of `pypolsys` project is to create an **easy to install and to deploy python package to start with Homotopy method**.

This wrapper has been developped to tackle polynomial systems obtained by high order perturbation of multi-parametric eigenvalue problems obtained with [`EasterEig`](https://github.com/nennigb/EasterEig). This yields to dense polynomial system. In this context, multivariate Horner algorithm has been added to quickly evaluate simultaneously the polynomial and its Jacobian matrix.

## POLSYS_PLP
`POLSYS_PLP` use a probability-one globally convergent homotopy method,
to find all finite **isolated complex solutions** to a system
F(X) = 0 of N polynomial equations in N unknowns with complex
coefficients.  A partitioned linear product (PLP) formulation is used
for the start system of the homotopy map.

The whole algorithm and its fortran90 implementation is described in
> Wise, Steven M., Andrew J. Sommese, and Layne T. Watson. "Algorithm 801: POLSYS_PLP: A partitioned linear product homotopy code for solving polynomial systems of equations." ACM Transactions on Mathematical Software (TOMS) 26.1 (2000): 176-200.

  - The sources are available on [netlib](http://www.netlib.org/opt/polsysplp.tgz).
  - The theoretical part of `POLSYS_PLP` is available [here](https://vtechworks.lib.vt.edu/bitstream/handle/10919/36933/thesis.pdf)
  - `POLSYS_PLP` is based on [`HOMPACK`](https://vtechworks.lib.vt.edu/bitstream/handle/10919/19612/TR-90-36.pdf)

To facilitate the build of this module, a copy of `POLSYS_PLP` *original* sources is included in the `pypolsys/801` folder.

## Installation

You'll need :
  * python (tested for v >= 3.5);
  * pip (optional);
  * fortran compiler (tested with `gfortran` and with `m2w64-toolchain` on windows)
  * lapack and blas installation, but the useful routines are also shipped with `POLSYS_PLP` sources. On debian based distribution it can be done with `sudo apt install libopenblas-dev`.

If needed, please see the steps given in the continuous integration script [workflows](.github/workflows/ci-ubuntu.yml).

### Using pip (preferred)
You can install `pypolsys` from pip:
```
pip install pypolsys [--user]
```
You can also install `pypolsys` after a download from github:
```
pip install path/to/pypolsys-version.tar.gz [--user]
```
or in _editable_ mode if you want to modify the sources
```
pip install -e path/to/pypolsys
```
after cloning the repos.

`pip` will install `numpy` (mandatory) and `sympy` and `scipy` (optional, required only for tests).

Note that on ubuntu, you will need to use `pip3` instead of `pip` and `python3` instead of `python`.

On linux, if you want to force the building with `POLSYS_PLP` lapack sources, you can tell it to pip with
```
pip install --no-use-pep517 --install-option="--buildLapack" path/to/pypolsys
```
The flag `--no-use-pep517` avoids side effect of `install-option` that disable wheels.
Using `POLSYS_PLP` lapack sources is the standard behavior on windows.

### Troubleshooting
With old `gfortran` version present on ubuntu 16.04 (see [here](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=84276)) the building may failed. To avoid this, `intent(in, out)` statement has to be removed from the original `pypolsys/801/polsys_plp.f90`. It seems to work out of the box with `gfortan-8`.

### Running tests

To execute the full test suite, run :
```
python -m pypolsys.test
```

## Usage
Note that this projet is a work in progress and the API may change.

### Get started
Let us consider the following example
```
x**2 + y + z - 1 = 0
x + y**2 + z - 1 = 0
x + y + z**2 - 1 = 0
```
> Cox, David, John Little, and Donal O'Shea. Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra. Springer Science & Business Media, 2013, From page 122.

With 5 distinct solutions in C : (1,0,0) (0,1,0), (0,0,1), (-1+√2, -1+√2, -1+√2), (-1 -√2, -1-√2, -1-√2). The main steps of `pypolsys`'s workflow are:
```python
import numpy as np
import pypolsys
# Declare the number of equation of the polynomial system
N = 3
# and the number of monoms of each equation (here 4 monoms for each equation).
n_coef_per_eq = np.array([4, 4, 4], dtype=np.int32)
# Provide the coefficients as 1D array
all_coef = np.array([1, 1,  1, -1,
                     1, 1,  1, -1,
                     1, 1,  1, -1], dtype=np.complex)
# then the degree of each monom
all_deg = np.zeros((np.sum(n_coef_per_eq), N), dtype=np.int32)
all_deg[0, 0] = 2
all_deg[1, 1] = 1
all_deg[2, 2] = 1
all_deg[4, 0] = 1
all_deg[5, 1] = 2
all_deg[6, 2] = 1
all_deg[8, 0] = 1
all_deg[9, 1] = 1
all_deg[10, 2] = 2
# Pass it to POLSYS_PLP
pypolsys.polsys.init_poly(N, n_coef_per_eq, all_coef, all_deg)
# Create homogeneous partition
# (N.B. Partitions are important to limit the number of paths
# to track associated to solutions at infinity)
part = pypolsys.utils.make_h_part(3)
# Pass it to POLSYS_PLP
pypolsys.polsys.init_partition(*part)
# Show coef
pypolsys.polsys.show_coef()
# Found 8 solutions, and track 8 paths; some solutions appear twice
bplp = pypolsys.polsys.solve(1e-8, 1e-14, 0.0)
# Get the roots, array of size (N+1) x bplp
r = pypolsys.polsys.myroots
# Get status of the solving process
pypolsys.polsys.report()
```
Note that due to projective transformation, the `r[N, :]` is the homogeneous variable. All results are returned ready to used (unscaled and untransformed). 

Another way, sometimes more convenient, is to setup the polynomial with function `pypolsys.utils.fromSympy`,

```python
import numpy as np
import pypolsys
import sympy as sym

x, y, z = sym.symbols('x, y, z')
pol = pypolsys.utils.fromSympy([sym.poly(x**2 + y + z - 1, (x, y, z)),
                                sym.poly(x + y**2 + z - 1, (x, y, z)),
                                sym.poly(x + y + z**2 - 1, (x, y, z))])
# Pass it to POLSYS_PLP
pypolsys.polsys.init_poly(*pol)
# Create homogeneous partition
part = pypolsys.utils.make_h_part(3)
# Pass it to POLSYS_PLP
pypolsys.polsys.init_partition(*part)
# Solve
bplp = pypolsys.polsys.solve(1e-8, 1e-14, 0.0)
# Get the roots, array of size (N+1) x bplp
r = pypolsys.polsys.myroots
# Get status of the solving process
pypolsys.polsys.report()
```
More examples are available in the [`pypolsys/test.py`](./pypolsys/test.py) file. For dense polynomial system, `dense` flag can be passed to the `solve` method to speed up evaluation.


### Organisation
This python package is divided into two parts
  - `pypolsys.polsys` corresponds to the interface to `POLSYS_PLP` and maps the following fortran subroutines (more doc is emmbedded in the f90 files) :
      * `init_poly`, initialize the polynomials monomials,
      * `init_partition`, initialize the variables partition,
      * `solve`, solve the polynomial system,
      * `refine`, refine the given paths,
      * `bezout`, compute the Bezout PLP number corresponding to the variable partition, usefull to test partition without solving,
      * `report`, show a report on the homotopy solving process
      * `show_partition` and `show_coef`, show the variables partitions and the monomials respectively.

      some attributs are also available :
      * `myroots`, containing the roots,
      * `path_status` and `solve_status`, containing the status of each tracked path and of the global solving process,
      * `num_jac_eval`, containing the number of evaluation of the Jacobian matrix.
  - `pypolsys.utils` contains utilities to facilitate the creation of the partition, polynomial or the estimation of the number of path (full python).

## Limitation
Due to shared variables in modules, only one polynomial instance can be handled. Concurrent calls to the `.so` file may create unexpected behavior.
To solve several polynomials, you can use `concurrent.futures`.

## How to contribute
This project started because we need some easy to deploy Homotopy solver. We are users but not experts on these methods.

If you want to contribute to `pypolsys`, your are welcomed! Don't hesitate to
  - report bugs, installation problems or ask questions on [issues](https://github.com/nennigb/pypolsys/issues);
  - propose some enhancements in the code or in documentation through **pull requests** (PR);
  - add or enhance support to other plateforms

To ensure code homogeneity among contributors, we use a source-code analyzer (eg. pylint).
Before submitting a PR, run the tests suite ;-)

## Developpement
If you need to modify `wrapper.f90` fortran source file. You will need to re-build the Fortran extension. This should be done into 2 steps
  1. The pyf file that avoid trouble with user_defined fortran type should be updated with `f2py` utilities :

    ```
    f2py -m polsys -h polsys.pyf wrapper.f90 --overwrite-signature
    ```
  2. Re run the install


## License
This file is part of pypolsys, a simple python wrapper to fortran package polsys_plp.
pypolsys is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
pypolsys is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with amc2moodle.  If not, see <https://www.gnu.org/licenses/>.
