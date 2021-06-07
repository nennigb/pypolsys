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

"""
Define some utility functions to pass *easily* polynomials and variables
partitions to `PPLSYS_PLP`.

Changing the partition has a strong impact on the number of paths to track
in the homotopy. `PPLSYS_PLP` supports
  * 1-homogeneous, denoted by h,
  * m-homogeneous, denoted by mh,
  * partitioned linear product, denoted by plp.

For more details on `PPLSYS_PLP` partition, See
> Wise, Steven M., Andrew J. Sommese, and Layne T. Watson. "Algorithm 801:
> POLSYS_PLP: A partitioned linear product homotopy code for solving polynomial
> systems of equations." ACM Transactions on Mathematical Software (TOMS)
> 26.1 (2000): 176-200.

For discussion on partitions and polynomials system types (dense, sparse, determinental), see
> Verschelde, Jan. "Polynomial homotopies for dense, sparse and determinantal
> systems." arXiv preprint math/9907060 (1999).
and the references therein.
"""

import numpy as np
from pypolsys import polsys

def fromSympy(P):
    """ Create polsys polynomial from sympy list of Poly. All variables
    should be present in the generator.

    Parameters
    ----------
    P : list
        A list of sympy poly objects.

    Returns
    -------
    N : int
        The number of variables and equations
    n_coef_per_eq : array
        The number of terms in each polynomial equation. Sum(n_coef_per_eq) give the
        total number of coefficeint
    all_coef : array
        Contains succesivelly all the coefficients of each monimial for each equation,
        starting by the first one.
    all_deg : array
        Contains succesivelly all the degrees of each monimial for each equation,
        starting by the first one.
        This is an array such the line number is the terms number and the column contain
        contains the degree of each variable.
        ex : if ther 4th term is (3.12+2j) x1^2 x2^3 x3^0, this yields
        all_deg(4, :) = [2, 3, 0] and all_coef(4) = 3.12+2j

    Remarks
    -------
    Use as a tuple, the output can be directly pass to `init_poly`.

    Examples
    --------
    >>> import sympy as sym
    >>> x, y = sym.symbols('x, y')
    >>> fromSympy([sym.poly(2*x + 3*y - x*y -3, (x,y)),\
                   sym.poly(3*x**2 + x*y - 1, (x,y))])  # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    (2,
    array([4, 3]...),
    array([-3.+0.j,  3.+0.j,  2.+0.j, -1.+0.j,
           -1.+0.j,  1.+0.j,  3.+0.j]),
    array([[0, 0],
           [0, 1],
           [1, 0],
           [1, 1],
           [0, 0],
           [1, 1],
           [2, 0]]...))
    """
    # Number of variables
    N = len(P)
    # Check if all polys have the same Generator
    if all([p.gens == P[0].gens for p in P]) is False:
        raise ValueError('All polynomial should have the same generator (.gens).')

    coef_list = []
    deg_list = []

    # Conversion
    n_coef_per_eq = np.zeros((N,), dtype=np.int32)
    for n, p in enumerate(P):
        pd = p.as_dict()
        k = 0
        coef_list_n = []
        deg_list_n = []
        for key, val in pd.items():
            coef_list_n.append(val)
            deg_list_n.append(list(key))
            k += 1
        # store the number of terms for this polynomial
        n_coef_per_eq[n] = k
        # convert to array
        coef_arr_n = np.array(coef_list_n, dtype=complex)
        deg_arr_n = np.array(deg_list_n, dtype=np.int32)
        # Sort Them
        index = np.lexsort(np.fliplr(deg_arr_n).T)
        deg_list.append(deg_arr_n[index, :].copy())
        coef_list.append(coef_arr_n[index].copy())

    return N, n_coef_per_eq, np.hstack(coef_list), np.vstack(deg_list)


def from1Darray(P):
    """ Create polsys polynomial from 1D array corresponding to a univariate
    polynomial. The first term is the constant; assume dense representation.

    Parameters
    ----------
    P : iterable
        Coefficient of the polynomial ordered by ascending order.

    Returns
    -------
    N : int
        The number of variables and equations
    n_coef_per_eq : array
        The number of terms in each polynomial equation. Sum(n_coef_per_eq) give the
        total number of coefficeint
    all_coef : array
        Contains succesivelly all the coefficients of each monimial for each equation,
        starting by the first one.
    all_deg : array
        Contains succesivelly all the degrees of each monimial for each equation,
        starting by the first one.
        This is an array such the line number is the terms number and the column contain
        contains the degree of each variable.
        ex : if ther 4th term is (3.12+2j) x1^2 x2^3 x3^0, this yields
        all_deg(4, :) = [2, 3, 0] and all_coef(4) = 3.12+2j

    Remarks
    -------
    Use as a tuple, the output can be directly pass to `init_poly`.

    Examples
    --------
    Consider the following example
        ```
        x**2 - 3*x + 2 = 0
        ```
    With 2 solutions in C : {2, 1}
    >>> from1Darray([2., -3., 1.])  # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    (1,
     array([3]...),
     array([ 2.+0.j, -3.+0.j,  1.+0.j]),
     array([[0],
            [1],
            [2]]...)
    """
    # Number of variables
    N = len(P)
    return (1, np.array([N], dtype=np.int32), np.array(P, dtype=complex),
            np.arange(0, N, dtype=np.int32).reshape(-1, 1))


def solve_univar(P, tracktol=1e-10, finaltol=1e-12, singtol=1e-14, dense=False):
    """ Solve univariate polynomial ordered by ascending power order.

    This function is a short hand to solve univariate case where the calling
    sequence can be simplified.

    Parameters
    ----------
    P : iterable
        Coefficient of the polynomial ordered by ascending order.
    tracktol : float
        is the local error tolerance allowed the path tracker along
        the path.
    finaltol : float
        is the accuracy desired for the final solution.  It is used
        for both the absolute and relative errors in a mixed error criterion.
   singtol : float
       is the singularity test threshold used by `SINGSYS_PLP`.  If
       `singtol <= 0.0` on input, then `singtol` is reset to a default value.
    dense : bool
       if `True`, select the `TARGET_SYSTEM_USER` optimized for dense polynomial
       (horner). If `False` (default) the defaut `POLSYS_PLP TARGET_SYSTEM` is
       used. The default choice is safer.

    Returns
    -------
    roots : array
        are the complex roots of the polynomial.

    Remarks
    --------
    For convenience the homogeneous variable due to the complex projective space
    are not return with this simplified interface.

    Examples
    --------
    Consider the following example
    ```
    x**2 - 3*x + 2 = 0
    ```
    With 2 solutions in C : {2, 1}
    >>> roots = solve_univar([2., -3., 1.])  # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    >>> np.sort(roots.real)
    array([1., 2.])
    """
    out4polsys = from1Darray(P)
    # Pass it to POLSYS_PLP
    polsys.init_poly(*out4polsys)
    # Create homogeneous partition
    part = make_h_part(1)
    # Pass it to POLSYS_PLP
    polsys.init_partition(*part)
    # Solve
    bplp = polsys.solve(tracktol, finaltol, singtol, dense)
    # Get the roots, array of size (N+1) x bplp
    roots = polsys.myroots
    return roots[0, :]



def di(all_deg, n_coef_per_eq):
    """ Return the degree of each equation and the total degree.

    This quantities are usefull to estimate the number of isolated solution
    of the polynomial system [Bezout Theorem].

    Parameters
    ----------
    all_deg : array
        Contains succesivelly all the degrees of each monimial for each equation,
        starting by the first one.

    n_coef_per_eq : array
        The number of terms in each polynomial equation. Sum(n_coef_per_eq) give the
        total number of coefficient

    Returns
    -------
    d : int
        The degree.
    dt : int
        the total degree.

    Examples
    --------
    >>> import sympy as sym
    >>> x, y = sym.symbols('x y')
    >>> N, n_coef_per_eq, all_coef, all_deg = fromSympy([sym.poly(2*x + 3*y - x*y -3, (x,y)),\
                                                         sym.poly(3*x**2 + x*y - 1, (x,y))])  # doctest: +NORMALIZE_WHITESPACE
    >>> di(all_deg, n_coef_per_eq)
    ([2, 2], 4)
    """
    d = []
    nini = 0
    for N in n_coef_per_eq:
        dmax = np.max(np.sum(all_deg[nini:nini+N, :], axis=1))
        d.append(dmax)
        nini += N

    dt = np.prod(d)

    return d, dt


def make_mh_part(N, var_list):
    """ Create `POLSYS_PLP` arguments for m-homogeneous partition.

    Multi-homogeneous or m-homogeneous partition allow to group unknwons into sets,
    The partition will be the same for all equations.

    Parameters
    ----------
    N : int
        The number of variables and equations.
    var_list : list
        A list containing the group of variables used for all equations.
        All the variables must be present only once in each sublist.
        e. g. with 3 variables: [[1], [2], [3]], or [[1, 2], [3]] or
        [[3, 2], [1]].
        Note the index **must start at 1** to be compatible with fortran.

    Returns
    -------
    N : int
        The number of variable and equations.
    num_set :  array(N)
       The number of set for each polynomial equation.
    num_indices : array(N,  NUM_SETS)
        Contains the 'list' of all the coefficients, sorting such poly1 coef, poly2 coef...
    index : array (N, NUM_SETS, NUM_INDICES)
        Contains the 'list' of all the variable, in each set for each equation.

    Remarks
    -------
    Use as a tuple, the output can be directly pass to `init_partition`.

    Examples
    --------
    >>> make_mh_part(3, [[1, 2], [3]]) # doctest: +NORMALIZE_WHITESPACE
    (3,
    array([2., 2., 2.]),
    array([[2, 1],
           [2, 1],
           [2, 1]]),
    array([[[1., 2., 0.],
            [3., 0., 0.]],
          [[1., 2., 0.],
           [3., 0., 0.]],
          [[1., 2., 0.],
           [3., 0., 0.]]]))
    """
    # multi-homogeneous order
    m = len(var_list)
    # populated matrix for fortran
    num_set = np.ones((N,)) * m
    num_indice = np.array([len(part) for part in var_list])
    num_indices = np.repeat(num_indice.reshape(1, -1), N, axis=0)
    # create full index array
    # this the same for all row
    mindex = np.zeros((1, m, N))
    for i, part in enumerate(var_list):
        mindex[0, i, 0:num_indice[i]] = np.array(part)
    index = np.repeat(mindex, N, axis=0)
    return N, num_set, num_indices, index


def make_h_part(N):
    """ Create `POLSYS_PLP` arguments for homogeneous partition.

    The partition will contain all variables, this equivalent to a 1-homogeneous
    partition. With this partition, the number of tracked paths will correspond to
    the total degree (Bézout Number).

    Parameters
    ----------
    N : int
        The number of variables and equations.

    Returns
    -------
    N : int
        The number of variable and equations.
    num_set :  array(N)
       The number of set for each polynomial equation.
    num_indices : array(N,  NUM_SETS)
        Contains the 'list' of all the coefficients, sorting such poly1 coef, poly2 coef...
    index : array (N, NUM_SETS, NUM_INDICES)
        Contains the 'list' of all the variable, in each set for each equation.

    Remarks
    -------
    Use as a tuple, the output can be directly pass to `init_partition`.
    """
    return make_mh_part(N, [list(range(1, N+1))])


def make_plp_part(N, var_list):
    """ Create `POLSYS_PLP` arguments for plp-homogeneous partition.

    Compared to m-homogeneous partition, plp partition allow to change
    the partition for each equation. With this partition, the number of
    tracked paths will correspond to Bézout PLP number.

    Parameters
    ----------
    N : int
        The number of variables and equations.
    var_list : list
        A list containing the group of variables used for each equation.
        All the variables must be present only once in each equation sublist.
        e. g. with 3 variables, len(var_list) = 3 and the sublists looks like
        [[[1], [2, 3]], [[2], [1, 3]], [[3], [1, 2]]] or
        [[[1, 2, 3]], [[2], [1], [3]], [[3], [1, 2]]]
        Note the index **must start at 1** to be compatible with fortran.

    Returns
    -------
    N : int
        The number of variable and equations.
    num_set :  array(N)
       The number of set for each polynomial equation.
    num_indices : array(N,  NUM_SETS)
        Contains the 'list' of all the coefficients, sorting such poly1 coef, poly2 coef...
    index : array (N, NUM_SETS, NUM_INDICES)
        Contains the 'list' of all the variable, in each set for each equation.

    Remarks
    -------
    Use as a tuple, the output can be directly pass to `init_partition`.
    """
    # Check if N == len(var_list)
    if N != len(var_list):
        raise ValueError('len(var_list) is {} instead of {}'.format(len(var_list), N))

    # populated input matrices for fortran
    num_set = np.array([len(part) for part in var_list])
    num_indices = np.zeros((N, N))
    # create full index array
    # this the same for all row
    index = np.zeros((N, max(num_set), N))
    for i, part in enumerate(var_list):
        for j, var in enumerate(part):
            num_indices[i, j] = len(var)
            index[i, j, 0:len(var)] = np.array(var)

    return N, num_set, num_indices, index


def toDense(N, n_coef_per_eq, all_coef, all_deg, preserve=True):
    """ Convert a *sparse* list of mononials suitable for `POLSYS_PLP`
    into a list of ordered monimials of a dense polynomials.

    Such preprocessing is **mandatory to use multiple variable Horner scheme**,
    ie `dense=True` in `polysys.solve`.
    For nearly dense polynomial, this method can significantly speed the polynomial
    and the Jacobian evaluation during the tracking.

      - Some absent monomial may be added to preserve the dense structure.
      - All monomials are listed left to right. The rightmost variable degree will
        change for each terms.
      - If a monimial is added to preserve the structure, its degree is set to 0,
        to avoid to change the Bezout Number computation in `POLSYS_PLP`.

    Parameters
    ----------
    N : int
        The number of variables and equations
    n_coef_per_eq : array
        The number of terms in each polynomial equation. Sum(n_coef_per_eq) give the
        total number of coefficeint
    all_coef : array
        Contains succesivelly all the coefficients of each monimial for each equation,
        starting by the first one.
    all_deg : array
        Contains succesivelly all the degrees of each monimial for each equation,
        starting by the first one.
        This is an array such the line number is the terms number and the column contain
        contains the degree of each variable.
        ex : if ther 4th term is (3.12+2j) x1^2 x2^3 x3^0, this yields
        all_deg(4, :) = [2, 3, 0] and all_coef(4) = 3.12+2j
    preserve : bool, optional
        Put the degree of vanishing coefficient to 0 to avoid to change Bezout number
        computation in `POLSYS_PLP` if dense structure is used. The default is True.

    Returns
    -------
    N : int
        The number of variables and equations
    n_coef_per_eq : array
        The number of terms in each polynomial equation. Sum(n_coef_per_eq) give the
        total number of coefficeint
    all_coef : array
        Contains succesivelly all the coefficients of each monimial for each equation,
        starting by the first one.
    all_deg : array
        Contains succesivelly all the degrees of each monimial for each equation,
        starting by the first one.
        This is an array such the line number is the terms number and the column contain
        contains the degree of each variable.

    Remarks
    -------
    Use as a tuple, the output can be directly pass to `init_poly`.
    Possible to enhance Horner speed if the highest degree is on the first variable.

    Examples
    --------
    >>> import sympy as sym
    >>> x, y = sym.symbols('x y')
    >>> sparse = fromSympy([sym.poly(3j*x**3 + 1, (x))])
    >>> toDense(*sparse)  # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
      > 200.0 % more coefs after conversion to dense.
    (1,
    array([4]...), array([1.+0.j, 0.+0.j, 0.+0.j, 0.+3.j]),
    array([[0],
          [0],
          [0],
         [3]]...))
    """
    # dmax contains the degree max in all variable for each equation
    dmax = np.zeros((N, N), dtype=np.int32)

    # setup new coefs matrix
    deg_list = []
    coef_list = []
    n_coef_per_eq_ = np.zeros_like(n_coef_per_eq)
    c = []
    # loop or equations
    start = 0
    for n in range(0, N):
        Nn = n_coef_per_eq[n]
        deg_n = all_deg[start:(start+Nn)]
        dmax[n, :] = np.max(deg_n, axis=0)
        # store ite
        deg_mat = np.zeros(dmax[n, :] + 1, dtype=complex)

        for coef, deg in zip(all_coef[start:(start+Nn)], deg_n):
            deg_mat[tuple(deg)] = coef
        c.append(deg_mat)
        # extract the indices and value
        deg_i = np.zeros((deg_mat.size, N), dtype=np.int32)
        coef_i = np.zeros((deg_mat.size,), dtype=complex)
        n_coef_per_eq_[n] = deg_mat.size
        k = 0
        for index, val in np.ndenumerate(deg_mat):
            if (val == 0) and preserve:
                # Put degree to 0 to avoid to change bezout Number for 0
                deg_i[k, :] = np.zeros_like(index)
            else:
                deg_i[k, :] = np.array(index)
            coef_i[k] = val
            k += 1
        deg_list.append(deg_i)
        coef_list.append(coef_i)

        start += Nn

    print('  > {} % more coefs after conversion to dense.'.format(np.sum(n_coef_per_eq_)
                                                         / np.sum(n_coef_per_eq)
                                                         * 100.))
    return N, n_coef_per_eq_, np.hstack(coef_list), np.vstack(deg_list)


# if __name__ == '__main__':
#     import doctest
#     doctest.testmod()
