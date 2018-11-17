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

from cext_common cimport _check_shape, _prepare_array, GBasis
from gbasis.pxds cimport c_gbasis

__all__ = [
    'GOBasis1',
]


#
# Internal business
#



cdef class GOBasis1(GBasis):
    # cdef public list _biblio
    # cdef c_gbasis.GOBasis* _this

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

    def compute_overlap(self, double[:, ::1] output=None) -> np.ndarray:
        r"""Compute the overlap integrals in a Gaussian orbital basis.

        .. math::
            \braket{\chi_{i}}{\chi_{j}}

        Parameters
        ----------
        output
            shape = (nbasis, nbasis). It will be overwritten.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.
        """
        output = _prepare_array(output, (self.nbasis, self.nbasis), 'output')
        self._this.compute_overlap(&output[0, 0])
        return np.asarray(output)

    def compute_kinetic(self, double[:, ::1] output=None) -> np.ndarray:
        r"""Compute the kinetic energy integrals in a Gaussian orbital basis.

        .. math::
            \mel{\chi_{i}}{\frac{\nabla^{2}}{2}}{\chi_{j}}

        Parameters
        ----------
        output
            shape = (nbasis, nbasis). It will be overwritten.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.
        """
        output = _prepare_array(output, (self.nbasis, self.nbasis), 'output')
        self._this.compute_kinetic(&output[0, 0])
        return np.asarray(output)

    def compute_nuclear_attraction(self, double[:, ::1] coordinates not None,
                                   double[::1] charges not None,
                                   double[:, ::1] output=None) -> np.ndarray:
        r"""Compute the nuclear attraction integral in a Gaussian orbital basis.

        .. math::
            \mel{\chi_{i}}{\frac{1}{\abs{r-R}}}{\chi_{j}}

        Parameters
        ----------
        coordinates
            A float array with shape (ncharge, 3) with Cartesian coordinates
            of point charges that define the external field.
        charges
            A float array with shape (ncharge,) with the values of the
            charges.
        output
            Two-index object, optional.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.
        """
        # type checking
        _check_shape(coordinates, (-1, 3), 'coordinates')
        natom = coordinates.shape[0]
        _check_shape(charges, (natom,), 'charges')
        output = _prepare_array(output, (self.nbasis, self.nbasis), 'output')
        # actual job
        self._this.compute_nuclear_attraction(
            &charges[0], &coordinates[0, 0], natom, &output[0, 0])
        return np.asarray(output)

    def compute_erf_attraction(self, double[:, ::1] coordinates not None,
                               double[::1] charges not None, double mu=0.0,
                               double[:, ::1] output=None) -> np.ndarray:
        r"""Compute the model nuclear attraction integral with the long-range potential

        .. math::
            \mel{\chi_{i}}{\frac{\mathrm{erf}(\mu r)}{r}}{\chi_{j}}

        Parameters
        ----------
        coordinates
            A float array with shape (ncharge, 3) with Cartesian coordinates
            of point charges that define the external field.
        charges
             A float array with shape (ncharge,) with the values of the
             charges.
        mu : float
             Parameter for the erf(mu r)/r potential. Default is zero.
        output
            Two-index object, optional.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.

        """
        # type checking
        _check_shape(coordinates, (-1, 3), 'coordinates')
        natom = coordinates.shape[0]
        _check_shape(charges, (natom,), 'charges')
        output = _prepare_array(output, (self.nbasis, self.nbasis), 'output')
        # actual job
        self._this.compute_erf_attraction(
            &charges[0], &coordinates[0, 0], natom, &output[0, 0], mu)
        return np.asarray(output)

    def compute_gauss_attraction(self, double[:, ::1] coordinates not None,
                                 double[::1] charges not None, double c=1.0,
                                 double alpha=1.0, double[:, ::1] output=None) -> np.ndarray:
        r"""Compute the model nuclear attraction with a Gaussian potential.

        .. math::
            \mel{\chi_{i}}{c\exp(-\alpha r^{2})}{\chi_{j}}

        Parameters
        ----------
        coordinates
            A float array with shape (ncharge, 3) with Cartesian coordinates
            of point charges that define the external field.
        charges
             A float array with shape (ncharge,) with the values of the
             charges.
        c : float
            Coefficient of the gaussian, default=1.0.
        alpha : float
            Exponential parameter of the gaussian, default=1.0.
        output
            Two-index object, optional.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.

        """
        # type checking
        _check_shape(coordinates, (-1, 3), 'coordinates')
        natom = coordinates.shape[0]
        _check_shape(charges, (natom,), 'charges')
        output = _prepare_array(output, (self.nbasis, self.nbasis), 'output')
        # actual job
        self._this.compute_gauss_attraction(
            &charges[0], &coordinates[0, 0], natom, &output[0, 0], c, alpha)
        return np.asarray(output)

    def compute_multipole_moment(self, long[::1] xyz, double[::1] center not None,
                                 double[:, ::1] output=None) -> np.ndarray:
        """Compute the (multipole) moment integrals in a Gaussian orbital basis.

        .. math::
            \mel{\chi_{i}}{(x-C_{x})^{l}(y-C_{y})^{m}(z-C_{z})^{n}}{\chi_{j}}

        Parameters
        ----------
        xyz : numpy-array of int, shape=(3,).
            A integer (long) numpy-array with shape (3,) with the powers of x,y,z in the
            integrals.
        center : np.ndarray, shape = (3,)
            A numpy array of shape (3,) with the center [C_x, C_y, C_z] around which the
            moment integral is computed.
        output
            Two-index object, optional.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.
        """
        # type checking
        _check_shape(center, (3,), 'center')
        _check_shape(xyz, (3,), 'xyz')
        if xyz[0] < 0 or xyz[1] < 0 or xyz[2] < 0 or xyz[0] + xyz[1] + xyz[2] <= 0:
            raise ValueError('Exponents of the Cartesian multipole operators must be '
                             'positive and may not sum to zero.')
        output = _prepare_array(output, (self.nbasis, self.nbasis), 'output')
        # actual job
        self._this.compute_multipole_moment(
            &xyz[0], &center[0], &output[0, 0])
        return np.asarray(output)


