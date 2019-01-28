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

from gbasis.cext_common cimport GBasis, _check_shape, _prepare_array
from gbasis.pxds cimport c_gbasis

__all__ = [
    'GOBasis2',
]

cdef class GOBasis2(GBasis):
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

    def compute_electron_repulsion(self, double[:, :, :, ::1] output=None) -> np.ndarray:
        r"""Compute electron-electron repulsion integrals.

        .. math::
            \mel{\chi_{i}\chi_{j}}{\frac{1}{\abs{r_{1}-r_{2}}}}{\chi_{k}\chi_{l}}

        Parameters
        ----------
        output
            A Four-index object, optional.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.

        """
        self._biblio.append(('valeev2014',
                    'the efficient implementation of four-center electron repulsion integrals'))
        output = _prepare_array(output, (self.nbasis, self.nbasis, self.nbasis, self.nbasis), 'output')
        self._this.compute_electron_repulsion(&output[0, 0, 0, 0])
        return np.asarray(output)

    def compute_erf_repulsion(self, double mu=0.0, double[:, :, :, ::1] output=None) -> np.ndarray:
        r"""Compute short-range electron repulsion integrals.

        .. math::
            \mel{\chi_{i}\chi_{j}}{\frac{1}{\abs{r_{1}-r_{2}}}}{\chi_{k}\chi_{l}}

        Parameters
        ----------
        mu : float
            Parameter for the erf(mu r)/r potential. Default is zero.
        output
            A Four-index object, optional.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.

        """
        self._biblio.append(('valeev2014',
                 'the efficient implementation of four-center electron repulsion integrals'))
        self._biblio.append(('ahlrichs2006',
                 'the methodology to implement various types of four-center integrals.'))
        output = _prepare_array(output, (self.nbasis, self.nbasis, self.nbasis, self.nbasis), 'output')
        self._this.compute_erf_repulsion(&output[0, 0, 0, 0], mu)
        return np.asarray(output)

    def compute_gauss_repulsion(self, double c=1.0, double alpha=1.0,
                                double[:, :, :, ::1] output=None) -> np.ndarray:
        r"""Compute gaussian repulsion four-center integrals.

        .. math::
            \mel{\chi_{i}\chi_{j}}{c\exp(-\alpha r^{2})}{\chi_{k}\chi_{l}}

        Parameters
        ----------
        c : float
            Coefficient of the gaussian.
        alpha : float
            Exponential parameter of the gaussian.
        output
            A Four-index object, optional.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.

        """
        self._biblio.append(('valeev2014',
                 'the efficient implementation of four-center electron repulsion integrals'))
        self._biblio.append(('ahlrichs2006',
                 'the methodology to implement various types of four-center integrals.'))
        self._biblio.append(('gill1996',
                 'four-center integrals with a Gaussian interaction potential.'))
        self._biblio.append(('toulouse2004',
                 'four-center integrals with a Gaussian interaction potential.'))
        output = _prepare_array(output, (self.nbasis, self.nbasis, self.nbasis, self.nbasis), 'output')
        self._this.compute_gauss_repulsion(&output[0, 0, 0, 0], c, alpha)
        return np.asarray(output)

    def compute_ralpha_repulsion(self, double alpha=-1.0,
                                 double[:, :, :, ::1] output=None) -> np.ndarray:
        r"""Compute r^alpha repulsion four-center integrals.

        The potential has the following form:

        .. math::
            \mel{\chi_{i}\chi_{j}}{r^{\alpha}}{\chi_{k}\chi_{l}}

        with :math:`\alpha > -3`.

        Parameters
        ----------
        alpha : float
            The power of r in the interaction potential.
        output
            A Four-index object, optional.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.


        """
        self._biblio.append(('valeev2014',
                 'the efficient implementation of four-center electron repulsion integrals'))
        self._biblio.append(('ahlrichs2006',
                 'the methodology to implement various types of four-center integrals.'))
        output = _prepare_array(output, (self.nbasis, self.nbasis, self.nbasis, self.nbasis), 'output')
        self._this.compute_ralpha_repulsion(&output[0, 0, 0, 0], alpha)
        return np.asarray(output)

    def compute_delta_repulsion(self, double[:, :, :, ::1] output=None,
                                      double[:, ::1] shift=None) -> np.ndarray:
        r"""Compute electron-electron repulsion integrals.

        .. math::
            \mel{\chi_{i}\chi_{j}}{\delta(\mathbf{r})}{\chi_{k}\chi_{l}}

        Parameters
        ----------
        shift
            A (4, 3) array of cartesian coordinates to shift each Gaussian center. Each row
            corresponds to a Gaussian center in physicist's notation.
        output
            A Four-index object, optional.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.

        """
        self.biblio.append(('valeev2014',
                    'the efficient implementation of four-center electron repulsion integrals'))
        output = _prepare_array(output, (self.nbasis, self.nbasis, self.nbasis, self.nbasis), 'output')
        if shift is not None:
            self._this.compute_delta_repulsion(&output[0, 0, 0, 0], &shift[0,0])
        else:
            self._this.compute_delta_repulsion(&output[0, 0, 0, 0])
        return np.asarray(output)

    def compute_intra_density(self, double[:, :, :, ::1] output=None,
                              double[:, ::1] point=None) -> np.ndarray:
        r"""Compute electron-electron repulsion integrals.

        .. math::
            \mel{\chi_{i}\chi_{j}}{\delta(\mathbf{r})}{\chi_{k}\chi_{l}}

        Parameters
        ----------
        output
            A Four-index object, optional.
        point
            The intracular coordinate.

        Returns
        -------
        np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis) dtype=float
            If output is provided it, it will be returned. Otherwise, a new array will be returned.

        """
        self.biblio.append(('valeev2014',
                    'the efficient implementation of four-center electron repulsion integrals'))
        _check_shape(point, (-1, 3), 'coordinates')
        if point is None:
            point = np.zeros((-1, 3))
        output = _prepare_array(output, (self.nbasis, self.nbasis, self.nbasis, self.nbasis), 'output')
        self._this.compute_intra_density(&output[0, 0, 0, 0], &point[0, 0])
        return np.asarray(output)
