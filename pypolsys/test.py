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

import unittest
import doctest
import pypolsys
import numpy as np
from scipy.optimize import linear_sum_assignment
import os
import sys


def sort_ref_found(ref, r):
    """ Sort reference solution and the found solution using
    linear sum assignment.

    Parameters
    ----------
    ref : array
        The reference solution.
    r : array
        The found solution.

    Returns
    -------
    err : float
        The total error.
    """
    # Create Cost Matrix for linear sum assignment between
    # This function can also solve a generalization of the classic assignment problem
    # where the cost matrix is rectangular. If it has more rows than columns,
    # then not every row needs to be assigned to a column, and vice versa.
    # all combinations
    Nroots = r.shape[1]
    Nref = ref.shape[1]
    Cost = np.zeros((Nroots, Nref), dtype=float)
    for index, c in np.ndenumerate(Cost):
        Cost[index] = np.linalg.norm(ref[:, index[1]] - r[0:-1, index[0]])
    row_ind, col_ind = linear_sum_assignment(Cost)
    err = Cost[row_ind, col_ind].sum()
    return err

class TestBasic(unittest.TestCase):
    """ Test suite
    """

    def test_POLSYS_PLP1(self):
        """ Test example #1 of original `POLSYS_PLP` test suite.

        Two quadrics, no solutions at infinity, two real solutions.'
        It corresponds to 801/OUTPUT.DAT results
        """
        # Tolerance to compare with ref. solution
        tol = 1e-8
        # Define the monoms
        N = 2
        n_coef_per_eq = np.array([6, 6], dtype=np.int32)
        all_coef = np.array([-9.80e-04, 9.78e+05, -9.80e+00, -2.35e+02, 8.89e+04, -1.00e+00,
                             -1.00e-02, -9.84e-01, -2.97e+01, 9.87e-03, -1.24e-01, -2.50e-01],
                            dtype=complex)
        all_deg = np.zeros((np.sum(n_coef_per_eq), N), dtype=np.int32)
        all_deg[0, 0] = 2
        all_deg[1, 1] = 2
        all_deg[2, 0] = 1; all_deg[2, 1] = 1
        all_deg[3, 0] = 1
        all_deg[4, 1] = 1
        all_deg[5, 1] = 0
        all_deg[6, 0] = 2
        all_deg[7, 1] = 2
        all_deg[8, 0] = 1; all_deg[8, 1] = 1
        all_deg[9, 0] = 1
        all_deg[10, 1] = 1

        # Pass Polynomial coefficients and degree to POLSYS_PLP
        pypolsys.polsys.init_poly(N, n_coef_per_eq, all_coef, all_deg)
        # Create homogeneous partition
        part = pypolsys.utils.make_h_part(2)
        # Pass it to POLSYS_PLP
        pypolsys.polsys.init_partition(*part)
        # Solve with POLSYS_PLP
        bplp = pypolsys.polsys.solve(1e-8, 1e-14, 0.0)
        # Get back the roots
        r = pypolsys.polsys.myroots
        # Check with the reference solution
        x1_ref = np.array([9.08921229615392E-02-1.11115036326264E-17j,  2.34233851959134E+03-5.78409592525803E-12j, 1.61478579234411E-02+1.68496955498882E+00j, 1.61478579234120E-02-1.68496955498883E+00j])
        x2_ref = np.array([-9.11497098197499E-02+2.31749408304374E-17j, -7.88344824094162E-01+1.94092735010767E-15j, 2.67994739614473E-04+4.42802993973665E-03j, 2.67994739614405E-04-4.42802993973668E-03j])
        self.assertEqual(bplp, 4)
        # print(np.linalg.norm(np.sort(x1_ref) - np.sort(r[0, :])) )
        self.assertTrue(np.linalg.norm(np.sort(x1_ref)
                                       - np.sort(r[0, :])) < tol)
        self.assertTrue(np.linalg.norm(np.sort(x2_ref)
                                       - np.sort(r[1, :])) < tol)

    def test_3vars_homogeneous_partition(self):
        """ Consider the following example
        ```
        x**2 + y + z - 1 = 0
        x + y**2 + z - 1 = 0
        x + y + z**2 - 1 = 0
        ```
        > Cox, David, John Little, and Donal OShea. Ideals, varieties, and algorithms: an
        > introduction to computational algebraic geometry and commutative algebra.
        > Springer Science & Business Media, 2013, From page 122.

        With 5 solutions in C : (1, 0, 0) (0, 1, 0), (0, 0,1 ), (-1+√2, -1+√2, -1+√2),
        (-1 -√2, -1-√2, -1-√2).
        """
        # Tolerance to compare with ref. solution
        tol = 1e-8
        # 3 equations, 3 variables
        N = 3
        # Each containing 4 monoms
        n_coef_per_eq = np.array([4, 4, 4], dtype=np.int32)
        # The coefficients are provided as 1D array
        all_coef = np.array([1, 1,  1, -1,
                             1, 1,  1, -1,
                             1, 1,  1, -1], dtype=complex)
        # Degree of each monom
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
        part = pypolsys.utils.make_h_part(3)
        # Pass it to POLSYS_PLP
        pypolsys.polsys.init_partition(*part)
        # Found 8 solution, and track 8 paths
        bplp = pypolsys.polsys.solve(1e-8, 1e-14, 0.0)
        # Get the roots, array of size (N+1) x bplp
        r = pypolsys.polsys.myroots
        # Create References solution.
        ref = np.zeros((3, 5))
        sol = [(1., 0., 0.),
               (0., 1., 0.),
               (0., 0, 1.),
               (-1. + np.sqrt(2.), -1. + np.sqrt(2.), -1. + np.sqrt(2.)),
               (-1. - np.sqrt(2.), -1. - np.sqrt(2.), -1. - np.sqrt(2.))]
        for i, s in enumerate(sol):
            ref[:, i] = s

        # Need to sort both solutions in the same way
        err = sort_ref_found(ref, r)
        self.assertTrue(err < tol)

    def test_Cox_refine(self):
        """ Use 1-homogeneous partition and refine computation.
        """
        import sympy as sym
        # Define tols
        tol = 1e-4
        tracktol = 1e-8
        finaltol = 1e-8
        singtol = 0.

        tracktolr = 1e-14
        finaltolr = 1e-14
        singtol = 0.

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
        # Solve with loose tolerances
        bplp = pypolsys.polsys.solve(tracktol, finaltol, finaltol)
        # Get the roots, array of size (N+1) x bplp
        r0 = pypolsys.polsys.myroots
        # Get status of the solving process
        pypolsys.polsys.report()
        # Create References solution.
        ref = np.zeros((3, 5))
        sol = [(1., 0., 0.),
               (0., 1., 0.),
               (0., 0, 1.),
               (-1. + np.sqrt(2.), -1. + np.sqrt(2.), -1. + np.sqrt(2.)),
               (-1. - np.sqrt(2.), -1. - np.sqrt(2.), -1. - np.sqrt(2.))]
        for i, s in enumerate(sol):
            ref[:, i] = s

        # Compute error with loose bounds
        err0 = sort_ref_found(ref, r0)
        # Refine one-half of the solutions
        path1 = np.array([1, 3, 5, 8])
        pypolsys.polsys.refine(path1, tracktolr, finaltolr, finaltol)
        path2 = np.array([2, 4, 6, 7])
        pypolsys.polsys.refine(path2, tracktolr, finaltolr, finaltol)
        # Get refine solution
        rr = pypolsys.polsys.myroots
        # Compute error with tight bounds
        errr = sort_ref_found(ref, rr)
        # Check error is lower than before
        self.assertTrue(errr*100 < err0)


    def test_univariate(self):
        """ Test univariate polynomial solve build from sympy.

        Consider the following example
        ```
        x**2 - 3*x + 2 = 0
        ```
        With 2 solutions in C : (2) (1)
        """
        # Tolerance to compare with ref. solution
        tol = 1e-8
        import sympy as sym
        x = sym.symbols('x')
        sparse = pypolsys.utils.fromSympy([sym.poly(x**2 - 3*x + 2, (x))])
        out4polsys = pypolsys.utils.toDense(*sparse)
        # Pass it to POLSYS_PLP
        pypolsys.polsys.init_poly(*out4polsys)
        # Create homogeneous partition
        part = pypolsys.utils.make_h_part(1)
        # Pass it to POLSYS_PLP
        pypolsys.polsys.init_partition(*part)
        # Solve
        bplp = pypolsys.polsys.solve(1e-8, 1e-14, 0.0, True)
        # Get the roots, array of size (N+1) x bplp
        r = pypolsys.polsys.myroots
        # Create References solution
        ref = np.array([1., 2.], dtype=complex)
        # Need to sort both solutions in the same way
        err = np.linalg.norm(ref - np.sort(r[0, :]))
        self.assertTrue(err < tol)
        self.assertEqual(bplp, 2)

    def test_horner_dense(self):
        """ Test horner, dense and standard approach on a simple test case,
        without solution at infnity.
        """

        import sympy as sym
        tol = 1e-4  # due to reference solution
        tracktol = 1e-8
        finaltol = 1e-14
        singtol = 0.
        # Example from https://nextjournal.com/saschatimme/homotopy-continuation-crash-course
        x, y = sym.symbols('x, y')
        monoms = pypolsys.utils.fromSympy([sym.poly(x**2 + 2*y, (x, y)),
                                           sym.poly(y**2 - 2, (x, y))])
        # 1-Homogeneous
        pypolsys.polsys.init_poly(*monoms)
        part = pypolsys.utils.make_h_part(2)
        pypolsys.polsys.init_partition(*part)
        # Solve and Store results
        b1 = pypolsys.polsys.solve(tracktol, finaltol, singtol)
        r1 = pypolsys.polsys.myroots.copy()

        # Use dense
        dmonoms = pypolsys.utils.toDense(*monoms)
        pypolsys.polsys.init_poly(*dmonoms)
        pypolsys.polsys.init_partition(*part)
        # Solve and Store results
        b1dense = pypolsys.polsys.solve(tracktol, finaltol, singtol, False)
        r1dense = pypolsys.polsys.myroots.copy()

        # Use dense AND Horner
        # Need to associate
        pypolsys.polsys.init_poly(*dmonoms)
        pypolsys.polsys.init_partition(*part)
        # Solve and Store results
        b1denseH = pypolsys.polsys.solve(tracktol, finaltol, singtol, True)
        r1denseH = pypolsys.polsys.myroots.copy()

        # Ref solution
        ref = np.array([[-2.51712e-32-1.68179j, 1.41421-3.69779e-32j],
                        [1.68179+1.54074e-33j, -1.41421-1.54074e-33j],
                        [2.51712e-32+1.68179j, 1.41421-3.69779e-32j],
                        [-1.68179-1.54074e-33j, -1.41421-1.54074e-33j]]).T
        # Need to sort both solutions in the same way
        err0 = sort_ref_found(ref, r1)        # standard
        err1 = sort_ref_found(ref, r1dense)   # dense
        err2 = sort_ref_found(ref, r1denseH)  # dense + horner
        self.assertTrue(err0 < tol)
        self.assertTrue(err1 < tol)
        self.assertTrue(err2 < tol)


class TestToyModel(unittest.TestCase):
    r""" More advanced Test case based on 3 dof toy model, with two parameters.

    The monomials and the reference solution are loaded from npz file.
        $$\begin{equation}        
        \left\lbrace
        \begin{aligned}
        - 2 \lambda - \mu - \nu + \left(\lambda + 2\right) \left(\lambda + \mu + 1\right) \left(\lambda + \nu + 1\right) - 2 &= 0, \\
        \left(\lambda + 2\right) \left(\lambda + \mu + 1\right) + \left(\lambda + 2\right) \left(\lambda + \nu + 1\right) + \left(\lambda + \mu + 1\right) \left(\lambda + \nu + 1\right) - 2 &= 0, \\
        6 \lambda + 2 \mu + 2 \nu + 8 &= 0.
        \end{aligned}
        \right.
        \end{equation}$$
    """

    @classmethod
    def setUpClass(cls):
        """ Load polynomial and reference solution from file.
        """
        # Load data
        input_file = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                    'examples/data/toy_model.npz'))
        data = np.load(input_file, allow_pickle=True)
        cls.out4polsys = data['out4polsys']
        cls.ref = data['ref'].T

    def test_1homogeneous(self):
        """ Use 1-homogeneous partition.
        """
        # Define tols
        tol = 1e-8
        tracktol = 1e-10
        finaltol = 1e-10
        singtol = 0.

        # 1-Homogeneous
        pypolsys.polsys.init_poly(*self.out4polsys)
        part = pypolsys.utils.make_h_part(3)
        pypolsys.polsys.init_partition(*part)
        # Solve and Store results
        b = pypolsys.polsys.solve(tracktol, finaltol, singtol)
        r = pypolsys.polsys.myroots
        # Need to sort both solutions in the same way
        err = sort_ref_found(self.ref, r)
        self.assertTrue(err < tol)

    def test_bezout(self):
        """ Use 1-homogeneous partition and compare with `bezout`.
        """
        # Define tols
        tol = 1e-8
        tracktol = 1e-10
        finaltol = 1e-10
        singtol = 0.

        # 1-Homogeneous
        pypolsys.polsys.init_poly(*self.out4polsys)
        part = pypolsys.utils.make_h_part(3)
        pypolsys.polsys.init_partition(*part)
        # Solve and compute Bezout
        b = pypolsys.polsys.solve(tracktol, finaltol, singtol)
        # Compute only Bezout
        B = pypolsys.polsys.bezout(singtol)
        # Check them
        self.assertEqual(b, B)

    def test_3_homogeneous(self):
        """ Use 3-homogeneous partition. Not optimal, but must find all solution.
        """
        # Define tols
        tol = 1e-8
        tracktol = 1e-10
        finaltol = 1e-10
        singtol = 0.

        # 1-Homogeneous
        pypolsys.polsys.init_poly(*self.out4polsys)
        part = pypolsys.utils.make_mh_part(3, [[1], [2], [3]])
        pypolsys.polsys.init_partition(*part)
        # Solve and Store results
        b = pypolsys.polsys.solve(tracktol, finaltol, singtol)
        r = pypolsys.polsys.myroots
        pypolsys.polsys.report()
        # Need to sort both solutions in the same way
        err = sort_ref_found(self.ref, r)
        self.assertTrue(err < tol)

    def test_plp_homogeneous(self):
        """ Use PLP partition. Not optimal, but must find all solution.
        """
        # Define tols
        tol = 1e-8
        tracktol = 1e-10
        finaltol = 1e-10
        singtol = 0.

        # PLP-partition
        pypolsys.polsys.init_poly(*self.out4polsys)
        part = pypolsys.utils.make_plp_part(3, [[[1], [2], [3]],
                                                [[1], [2, 3]],
                                                [[1, 2, 3]]])
        pypolsys.polsys.init_partition(*part)
        # Solve and Store results
        b = pypolsys.polsys.solve(tracktol, finaltol, singtol)
        r = pypolsys.polsys.myroots
        # Need to sort both solutions in the same way
        err = sort_ref_found(self.ref, r)
        self.assertTrue(err < tol)

    def test_3_homogeneousH(self):
        """ Use 3-homogeneous partition. Not optimal, but must find all solution.
        """
        # Define tols
        tol = 1e-8
        tracktol = 1e-10
        finaltol = 1e-10
        singtol = 0.

        # 1-Homogeneous
        dense = pypolsys.utils.toDense(*self.out4polsys)
        pypolsys.polsys.init_poly(*dense)
        part = pypolsys.utils.make_mh_part(3, [[1], [2], [3]])
        pypolsys.polsys.init_partition(*part)
        # Solve and Store results
        b = pypolsys.polsys.solve(tracktol, finaltol, singtol, True)
        r = pypolsys.polsys.myroots
        # Need to sort both solutions in the same way
        err = sort_ref_found(self.ref, r)
        self.assertTrue(err < tol)


if __name__ == '__main__':
    # Run unittest test suite
    print('> Running tests...')
    loadTestsFromTestCase = unittest.defaultTestLoader.loadTestsFromTestCase
    suite = loadTestsFromTestCase(TestBasic)
    suite.addTest(loadTestsFromTestCase(TestToyModel))
    suite.addTest(doctest.DocTestSuite(pypolsys.utils))
    # define the runner
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    # runner doesn't change exit status
    if result.wasSuccessful():
        sys.exit(0)
    else:
        sys.exit(1)
