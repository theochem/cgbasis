# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --

#!python
#cython: embedsignature=True, language_level=3

"""C++ extensions"""
cimport numpy as np
import numpy as np

np.import_array()

from gbasis.cext_common cimport GBasis, _check_shape, _prepare_array, _get_shell_nbasis

from gbasis.pxds cimport c_gbasis
from gbasis.pxds.grids cimport c_fns
from gbasis.pxds.grids cimport c_nucpot

__all__ = [
    # common
    'GOBasisGrid',
    # nucpot
    'compute_grid_nucpot',

    # fns
    '_GB1DMGridDensityFn', '_GB1DMGridGradientFn', '_GB1DMGridGGAFn',
    '_GB1DMGridKineticFn', '_GB1DMGridHessianFn', '_GB1DMGridMGGAFn',
]



cdef class GOBasisGrid(GBasis):
    cdef public list _biblio
    cdef c_gbasis.GOBasis* _this

    def __cinit__(self, centers, shell_map, nprims, shell_types, alphas, con_coeffs):
        self._this = new c_gbasis.GOBasis(
            <double*>self._centers.data, <long*>self._shell_map.data,
            <long*>self._nprims.data, <long*>self._shell_types.data,
            <double*>self._alphas.data, <double*>self._con_coeffs.data,
            self._centers.shape[0], self._shell_types.shape[0], self._alphas.shape[0])
        self._baseptr = <c_gbasis.GBasis*> self._this

    def check_coeffs(self, double[:, ::1] coeffs, nocc=None):
        if coeffs.shape[0] != self.nbasis:
            raise TypeError('coeffs.shape[0] should be equal to nbasis. '
                            'Got {}. Expecting {}.'.format(coeffs.shape[0], self.nbasis))
        if coeffs.shape[1] > self.nbasis:
            raise TypeError('coeffs.shape[1] should be less than equal to nbasis. '
                            'Got {}. Expecting <={}.'.format(coeffs.shape[0], self.nbasis))
        if nocc is not None and coeffs.shape[1] < nocc:
            raise TypeError('coeffs.shape[1] should not be below the number of occupation '
                            'numbers. Got {}. Expected at least {}.'.format(coeffs.shape[1], nocc))


    # low-level compute routines, for debugging only
    def _compute_grid_point1(self, double[::1] output not None,
                            double[::1] point not None,
                            _GB1DMGridFn grid_fn not None):
        assert output.shape[0] == self.nbasis
        assert point.shape[0] == 3
        self._baseptr.compute_grid_point1(&output[0], &point[0], grid_fn._baseptr)

    def compute_grid_orbitals_exp(self, double[:, ::1] coeffs, double[:, ::1] points not None,
                                  long[::1] iorbs not None,
                                  double[:, ::1] output=None):
        r"""Compute the orbitals on a grid for a given set of expansion coefficients.

        **Warning:** the results are added to the output array!

        Parameters
        ----------
        coeffs : np.ndarray, shape=(nbasis, nfn), dtype=float
            The orbitals coefficients
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        iorbs : np.ndarray, shape=(n,), dtype=int
            Indexes of the orbitals to be computed. When not given, the orbitals with a
            non-zero occupation number are computed.
        output : np.ndarray, shape=(npoint, n), dtype=float
            An output array. The results are added to this array. When not given, an
            output array is allocated.

        Returns
        -------
        np.ndarray, shape=(npoint, n), dtype=float
            the output array. (It is allocated when not given.)
        """
        # Do some type checking
        self.check_coeffs(coeffs)
        nfn = coeffs.shape[1]
        _check_shape(points, (-1, 3), 'points')
        npoint = points.shape[0]
        norb = iorbs.shape[0]
        output = _prepare_array(output, (npoint, norb), 'output')
        # compute
        self._this.compute_grid1_exp(
            nfn, &coeffs[0, 0], npoint, &points[0, 0],
            norb, &iorbs[0], &output[0, 0])
        return np.asarray(output)

    def compute_grid_orb_gradient_exp(self, double[:, ::1] coeffs, double[:, ::1] points not None,
                                      long[::1] iorbs not None,
                                      double[:, :, ::1] output=None):
        r"""Compute the orbital gradient on a grid for a given set of expansion coefficients.

        **Warning:** the results are added to the output array!

        Parameters
        ----------
        coeffs : np.ndarray, shape=(nbasis, nfn), dtype=float
            Orbital coefficients
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        iorbs : np.ndarray, shape=(n,), dtype=int
            Indexes of the orbitals to be computed.
        output : np.ndarray, shape=(npoint, n, 3), dtype=float
            An output array. The results are added to this array. When not given, an
            output array is allocated.

        Returns
        -------
        output : np.ndarray, shape=(npoint, n, 3), dtype=float
            the output array. (It is allocated when not given.)
        """
        # Do some type checking
        self.check_coeffs(coeffs)
        nfn = coeffs.shape[1]
        _check_shape(points, (-1, 3), 'points')
        npoint = points.shape[0]
        norb = iorbs.shape[0]
        output = _prepare_array(output, (npoint, norb, 3), 'output')
        # compute
        self._this.compute_grid1_grad_exp(
            nfn, &coeffs[0, 0], npoint, &points[0, 0],
            norb, &iorbs[0], &output[0, 0, 0])
        return np.asarray(output)

    def _compute_grid1_dm(self, double[:, ::1] dm not None, double[:, ::1] points not None,
                          _GB1DMGridFn grid_fn not None, double[:, ::1] output not None,
                          double epsilon=0):
        """Compute some density function on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        grid_fn : _GB1DMGridFn
            Implements the function to be evaluated on the grid.
        output : np.ndarray, shape=(npoint, n), dtype=float
            Output array. The second dimension depends on grid_fn.
        epsilon : float
            Allow errors on the density of this magnitude for the sake of
            efficiency. Some grid_fn implementations may ignore this.

        Returns
        -------
        np.ndarray, shape=(npoint, n) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.
        """
        # Check the array shapes
        _check_shape(dm, (self.nbasis, self.nbasis,), 'dm')
        _check_shape(points, (-1, 3), 'points')
        npoint = points.shape[0]
        _check_shape(output, (npoint, grid_fn.dim_output), 'output')
        # Get the maximum of the absolute value over the rows
        cdef double[:] dmmaxrow = np.abs(dm).max(axis=0)
        # Go!
        self._this.compute_grid1_dm(
            &dm[0, 0], npoint, &points[0, 0], grid_fn._baseptr, &output[0, 0], epsilon,
            &dmmaxrow[0])

    def compute_grid_density_dm(self, double[:, ::1] dm not None,
                                double[:, ::1] points not None, double[::1] output=None,
                                double epsilon=0):
        """Compute the electron density on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint,), dtype=float
            Output array. When not given, it is allocated and returned.
        epsilon : float
            Allow errors on the density of this magnitude for the sake of
            efficiency. Some grid_fn implementations may ignore this.

        Returns
        -------
        output : np.ndarray, shape=(npoint,), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros(points.shape[0])
        self._compute_grid1_dm(dm, points, _GB1DMGridDensityFn(self.max_shell_type),
                               output[:, None], epsilon)
        return np.asarray(output)

    def compute_grid_gradient_dm(self, double[:, ::1] dm not None,
                                 double[:, ::1] points not None,
                                 double[:, ::1] output=None):
        """Compute the electron density gradient on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint, 3), dtype=float
            Output array. When not given, it is allocated and returned.

        Returns
        -------
        output : np.ndarray, shape=(npoint, 3), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0], 3), float)
        self._compute_grid1_dm(dm, points, _GB1DMGridGradientFn(self.max_shell_type), output)
        return np.asarray(output)

    def compute_grid_gga_dm(self, double[:, ::1] dm not None,
                            double[:, ::1] points not None,
                            double[:, ::1] output=None):
        """Compute the electron density and gradient on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint, 4), dtype=float
            Output array. When not given, it is allocated and returned. The first column
            contains the density. The last three columns contain the gradient.

        Returns
        -------
        output : np.ndarray, shape=(npoint, 4), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0], 4), float)
        self._compute_grid1_dm(dm, points, _GB1DMGridGGAFn(self.max_shell_type), output)
        return np.asarray(output)

    def compute_grid_kinetic_dm(self, double[:, ::1] dm not None,
                                double[:, ::1] points not None,
                                double[::1] output=None):
        """Compute the positive definite kinetic energy density on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint,), dtype=float
            Output array. When not given, it is allocated and returned.

        Returns
        -------
        output : np.ndarray, shape=(npoint,), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0],), float)
        self._compute_grid1_dm(dm, points, _GB1DMGridKineticFn(self.max_shell_type), output[:, None])
        return np.asarray(output)

    def compute_grid_hessian_dm(self, double[:, ::1] dm not None,
                                double[:, ::1] points not None,
                                double[:, ::1] output=None):
        """Compute the electron density Hessian on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint, 6), dtype=float
            Output array. When not given, it is allocated and returned. The columns are
            assigned as follows:

            * 0: element (0, 0) of the Hessian
            * 1: element (0, 1) of the Hessian
            * 2: element (0, 2) of the Hessian
            * 3: element (1, 1) of the Hessian
            * 4: element (1, 2) of the Hessian
            * 5: element (2, 2) of the Hessian

        Returns
        -------
        output : np.ndarray, shape=(npoint, 6), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0], 6), float)
        self._compute_grid1_dm(dm, points, _GB1DMGridHessianFn(self.max_shell_type), output)
        return np.asarray(output)

    def compute_grid_mgga_dm(self, double[:, ::1] dm not None,
                             double[:, ::1] points not None,
                             double[:, ::1] output=None):
        """Compute the MGGA quantities for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        This includes the density, the gradient, the Laplacian and the kinetic energy
        density.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint, 6), dtype=float
            Output array. When not given, it is allocated and returned. The assignment of
            the columns is as follows:

            * 0: density
            * 1: gradient x
            * 2: gradient y
            * 3: gradient z
            * 4: laplacian
            * 5: kinetic energy density

        Returns
        -------
        output : np.ndarray, shape=(npoint, 6), dtype=float
            The output array.
        """
        if output is None:
            output = np.zeros((points.shape[0], 6), float)
        self._compute_grid1_dm(dm, points, _GB1DMGridMGGAFn(self.max_shell_type), output)
        return np.asarray(output)

    def compute_grid_hartree_dm(self, double[:, ::1] dm not None,
                                double[:, ::1] points not None,
                                double[::1] output=None):
        """Compute the Hartree potential on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Density matrix, assumed to be symmetric.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint,), dtype=float
            Output array. When not given, it is allocated and returned.

        Returns
        -------
        output : np.ndarray, shape=(npoint,), dtype=float
            The output array.
        """
        # type checking
        _check_shape(dm, (self.nbasis, self.nbasis), 'dm')
        _check_shape(points, (-1, 3), 'points')
        npoint = points.shape[0]
        output = _prepare_array(output, (npoint,), 'output')
        # compute
        self._this.compute_grid2_dm(
            &dm[0, 0], npoint, &points[0, 0], &output[0])
        return np.asarray(output)

    def compute_grid_esp_dm(self, double[:, ::1] dm not None,
                            double[:, ::1] coordinates not None, double[::1] charges not None,
                            double[:, ::1] points not None,
                            double[::1] output=None):
        """Compute the electrostatic potential on a grid for a given density matrix.

        **Warning:** the results are added to the output array! This may be useful to
        combine results from different spin components.

        Parameters
        ----------
        dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Density matrix, assumed to be symmetric.
        coordinates : np.ndarray, shape=(natom, 3), dtype=float
            Cartesian coordinates of the atoms.
        charges : np.ndarray, shape=(natom,), dtype=float
            Atomic charges.
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        output : np.ndarray, shape=(npoint,), dtype=float
            Output array. When not given, it is allocated and returned.

        Returns
        -------
        output : np.ndarray, shape=(npoint,), dtype=float
            The output array.
        """
        output = self.compute_grid_hartree_dm(dm, points, output)
        cdef np.ndarray[ndim=1, dtype=double] tmp = np.asarray(output)
        tmp *= -1
        compute_grid_nucpot(coordinates, charges, points, output)
        return tmp

    def _compute_grid1_fock(self, double[:, ::1] points not None, double[::1] weights not None,
                            double[:, :] pots not None, _GB1DMGridFn grid_fn not None,
                            double[:, ::1] fock=None):
        """Compute a Fock operator from a some sort of potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, n)
            Derivative of the energy toward the density-related quantities
            at all grid points. The number of columns depends on grid_fn.
        grid_fn : _GB1DMGridFn
            Implements the function to be evaluated on the grid.
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Output two-index object.

        Returns
        -------
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            If fock is provided, it will be **added** to. If not, a new array is allocated and
            returned.
        """
        fock = _prepare_array(fock, (self.nbasis, self.nbasis), 'fock')
        _check_shape(points, (-1, 3), 'points')
        npoint = points.shape[0]
        _check_shape(weights, (npoint,), 'weights')
        _check_shape(pots, (npoint, grid_fn.dim_output), 'pots')
        if pots.strides[0] % 8 != 0:
            raise TypeError('stride[0] of the pots argument must be a multiple of 8.')
        if pots.strides[1] % 8 != 0:
            raise TypeError('stride[1] of the pots argument must be a multiple of 8.')
        pot_stride = (pots.strides[0]//8)
        if pots.shape[1] > 1:
            pot_stride *= (pots.strides[1]//8)
        self._this.compute_grid1_fock(
            npoint, &points[0, 0], &weights[0], pot_stride, &pots[0, 0], grid_fn._baseptr,
            &fock[0, 0])
        return fock

    def compute_grid_density_fock(self, double[:, ::1] points not None,
                                  double[::1] weights not None, double[:] pots not None,
                                  double[:, ::1] fock=None):
        """Compute a Fock operator from a density potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint,), dtype=float
            Derivative of the energy toward the density at all grid points.
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Output two-index object, optional.

        Returns
        -------
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            If fock is provided, it will be **added** to. If not, a new array is allocated and
            returned.
        """
        fock = self._compute_grid1_fock(
            points, weights, pots[:, None], _GB1DMGridDensityFn(self.max_shell_type), fock)
        return np.asarray(fock)

    def compute_grid_gradient_fock(self, double[:, ::1] points not None,
                                   double[::1] weights not None, double[:, :] pots not None,
                                   double[:, ::1] fock=None):
        """Compute a Fock operator from a density gradient potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, 3), dtype=float
            Derivative of the energy toward the density gradient components at all grid
            points.
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Output two-index object, optional.

        Returns
        -------
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            If fock is provided, it will be **added** to. If not, a new array is allocated and
            returned.
        """
        fock = self._compute_grid1_fock(
            points, weights, pots, _GB1DMGridGradientFn(self.max_shell_type), fock)
        return np.asarray(fock)

    def compute_grid_gga_fock(self, double[:, ::1] points not None,
                              double[::1] weights not None, double[:, :] pots not None,
                              double[:, ::1] fock=None):
        """Compute a Fock operator from GGA potential data.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, 4), dtype=float
            Derivative of the energy toward GGA ingredients (density and gradient) at all
            grid points.
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Output two-index object, optional.

        Returns
        -------
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            If fock is provided, it will be **added** to. If not, a new array is allocated and
            returned.
        """
        fock = self._compute_grid1_fock(
            points, weights, pots, _GB1DMGridGGAFn(self.max_shell_type), fock)
        return np.asarray(fock)

    def compute_grid_kinetic_fock(self, double[:, ::1] points not None,
                                  double[::1] weights not None, double[:] pots not None,
                                  double[:, ::1] fock=None) :
        """Compute a Fock operator from a kinetic-energy-density potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint,), dtype=float
            Derivative of the energy toward the kinetic energy density at all grid points.
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Output two-index object, optional.

        Returns
        -------
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            If fock is provided, it will be **added** to. If not, a new array is allocated and
            returned.
        """
        fock = self._compute_grid1_fock(
            points, weights, pots[:, None], _GB1DMGridKineticFn(self.max_shell_type), fock)
        return np.asarray(fock)

    def compute_grid_hessian_fock(self, double[:, ::1] points not None,
                                  double[::1] weights not None, double[:, :] pots not None,
                                  double[:, ::1] fock=None):
        """Compute a Fock operator from a density hessian potential.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, 6), dtype=float
            Derivative of energy toward components of the Hessian at all grid points. The
            columns are assigned as follows:

            * 0: element (0, 0) of the Hessian
            * 1: element (0, 1) of the Hessian
            * 2: element (0, 2) of the Hessian
            * 3: element (1, 1) of the Hessian
            * 4: element (1, 2) of the Hessian
            * 5: element (2, 2) of the Hessian

        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Output two-index object, optional.

        Returns
        -------
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            If fock is provided, it will be **added** to. If not, a new array is allocated and
            returned.
        """
        fock = self._compute_grid1_fock(
            points, weights, pots, _GB1DMGridHessianFn(self.max_shell_type), fock)
        return np.asarray(fock)

    def compute_grid_mgga_fock(self, double[:, ::1] points not None,
                               double[::1] weights not None, double[:, :] pots not None,
                               double[:, ::1] fock=None):
        """Compute a Fock operator from MGGA potential data.

        **Warning:** the results are added to the Fock operator!

        Parameters
        ----------
        points : np.ndarray, shape=(npoint, 3), dtype=float
            Cartesian grid points.
        weights : np.ndarray, shape=(npoint,), dtype=float
            Integration weights.
        pots : np.ndarray, shape=(npoint, 6), dtype=float
            Derivative of the energy toward density, gradient, Laplacian and kinetic
            energy density. The assignment of the columns is as follows:

            * 0: density
            * 1: gradient x
            * 2: gradient y
            * 3: gradient z
            * 4: laplacian
            * 5: kinetic energy density

        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            Output two-index object, optional.

        Returns
        -------
        fock : np.ndarray, shape=(nbasis, nbasis), dtype=float
            If fock is provided, it will be **added** to. If not, a new array is allocated and
            returned.
        """
        fock = self._compute_grid1_fock(
            points, weights, pots, _GB1DMGridMGGAFn(self.max_shell_type), fock)
        return np.asarray(fock)


#
# nucpot.cpp
#


def compute_grid_nucpot(double[:, ::1] coordinates not None,
                        double[::1] charges not None,
                        double[:, ::1] points not None,
                        double[::1] output not None) -> np.ndarray:
    """Compute the potential due to a set of (nuclear) point charges

    Parameters
    ----------
    coordinates
        A (N, 3) float numpy array with Cartesian coordinates of the
        atoms.
    charges
        A (N,) numpy vector with the atomic charges.
    points
        An (M, 3) array with grid points where the potential must be
        computed.
    output
        An (M,) output array in which the potential is stored.

    Returns
    -------
    output : np.ndarray, shape=(M, ), dtype=float
    """
    # type checking
    assert coordinates.shape[1] == 3
    ncharge = coordinates.shape[0]
    assert charges.shape[0] == ncharge
    assert points.shape[1] == 3
    npoint = points.shape[0]
    assert output.shape[0] == npoint
    # actual computation
    c_nucpot.compute_grid_nucpot(
        &coordinates[0,0], &charges[0], ncharge,
        &points[0,0], &output[0], npoint)



#
# fns wrappers (for testing and use in this module)
#


cdef class _GB1DMGridFn:
    """Wrapper for fns.GB1DMGridFn, for testing only"""
    cdef c_fns.GB1DMGridFn* _baseptr

    def __dealloc__(self):
        del self._baseptr

    @property
    def nwork(self):
        return self._baseptr.get_nwork()

    @property
    def max_shell_type(self):
        return self._baseptr.get_max_shell_type()

    @property
    def max_nbasis(self):
        return self._baseptr.get_max_nbasis()

    @property
    def shell_type0(self):
        return self._baseptr.get_shell_type0()

    @property
    def dim_work(self):
        return self._baseptr.get_dim_work()

    @property
    def dim_output(self):
        return self._baseptr.get_dim_output()

    def reset(self, long shell_type0, double[::1] r0 not None, double[::1] point not None):
        assert r0.shape[0] == 3
        assert point.shape[0] == 3
        self._baseptr.reset(shell_type0, &r0[0], &point[0])

    def add(self, double coeff, double alpha0,
            double[::1] scales0 not None):
        assert scales0.shape[0] == _get_shell_nbasis(abs(self._baseptr.get_shell_type0()))
        self._baseptr.add(coeff, alpha0, &scales0[0])

    def cart_to_pure(self):
        self._baseptr.cart_to_pure()

    def get_work(self, shape0):
        """This returns a **copy** of the c++ work array.

           Returning a numpy array with a buffer created in c++ is dangerous.
           If the c++ array becomes deallocated, the numpy array may still
           point to the deallocated memory. For that reason, a copy is returned.
           Speed is not an issue as this class is only used for testing.
        """
        cdef np.npy_intp shape[2]
        assert shape0 > 0
        assert shape0 <= self.max_nbasis
        shape[0] = shape0
        if self.dim_work == 1:
            tmp = np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, <void*> self._baseptr.get_work())
        else:
            shape[1] = self.dim_work
            tmp = np.PyArray_SimpleNewFromData(2, shape, np.NPY_DOUBLE, <void*> self._baseptr.get_work())
        return tmp.copy()


cdef class _GB1DMGridDensityFn(_GB1DMGridFn):
    cdef c_fns.GB1DMGridDensityFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new c_fns.GB1DMGridDensityFn(max_nbasis)
        self._baseptr = <c_fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridGradientFn(_GB1DMGridFn):
    cdef c_fns.GB1DMGridGradientFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new c_fns.GB1DMGridGradientFn(max_nbasis)
        self._baseptr = <c_fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridGGAFn(_GB1DMGridFn):
    cdef c_fns.GB1DMGridGGAFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new c_fns.GB1DMGridGGAFn(max_nbasis)
        self._baseptr = <c_fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridKineticFn(_GB1DMGridFn):
    cdef c_fns.GB1DMGridKineticFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new c_fns.GB1DMGridKineticFn(max_nbasis)
        self._baseptr = <c_fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridHessianFn(_GB1DMGridFn):
    cdef c_fns.GB1DMGridHessianFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new c_fns.GB1DMGridHessianFn(max_nbasis)
        self._baseptr = <c_fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridMGGAFn(_GB1DMGridFn):
    cdef c_fns.GB1DMGridMGGAFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new c_fns.GB1DMGridMGGAFn(max_nbasis)
        self._baseptr = <c_fns.GB1DMGridFn*> self._this
