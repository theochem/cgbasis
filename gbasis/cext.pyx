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
#cython: embedsignature=True

"""C++ extensions"""
from typing import Iterable, Tuple, Type, List, Union

cimport numpy as np
import numpy as np

np.import_array()

from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

from numpy_wrapper cimport PyArray_ENABLEFLAGS

cimport common
cimport gbasis
cimport ints
cimport fns


cimport cholesky
cimport gbw
cimport nucpot

import atexit

__all__ = [
    # common
    '_get_shell_nbasis',
    # gbasis
    '_gob_cart_normalization', '_gob_pure_normalization',
    'GOBasis',
    # ints
    '_GB4Integral',
    '_GB4ElectronRepulsionIntegralLibInt',
    '_GB4ErfIntegralLibInt', '_GB4GaussIntegralLibInt',
    '_GB4RAlphaIntegralLibInt',
    '_GB4DeltaIntegralLibInt',
    '_GB4IntraDensIntegralLibInt',

    # fns
    '_GB1DMGridDensityFn', '_GB1DMGridGradientFn', '_GB1DMGridGGAFn',
    '_GB1DMGridKineticFn', '_GB1DMGridHessianFn', '_GB1DMGridMGGAFn',


]


#
# Internal business
#

cdef _check_shape(array: Union[np.ndarray, memoryview], shape: Tuple[int], name: str):
    """Checks the shape of an array.

    Parameters
    ----------
    array
        Array to be checked
    shape
        Expected shape. Negative sizes are not checked.
    name
        The name of the array, used in the error messages.

    Raises
    ------
    TypeError
        When the dimension or the shape of the array differs from the expected shape.
    """
    if array.ndim != len(shape):
        raise TypeError('Array \'{}\' has ndim={}. Expecting ndim={}.'.format(
            name, array.ndim, len(shape)))
    for i, si in enumerate(shape):
        if 0 <= si != array.shape[i]:
            raise TypeError('Array \'{}\' has shape[{}]={}. Expecting shape[{}]={}.'.format(
            name, i, array.shape[i], i, si))


cdef _prepare_array(array: Union[np.ndarray, None], shape: Tuple[int], name: str):
    """Check array shape or initialize it if None.

    Parameters
    ----------
    array
        Array to be checked (or initialized if None).
    shape
        The expected shape.
    name
        The name of the array, used in the error messages.

    Returns
    -------
    np.ndarray

    Raises
    ------
    TypeError
        When the dimension or the shape of the array differs from the expected shape.
    """
    if array is None:
        array = np.zeros(shape)
    else:
        _check_shape(array, shape, name)
    return array

#
# common.cpp wrappers
#


cpdef _get_shell_nbasis(long shell_type):
    result = common.get_shell_nbasis(shell_type)
    if result <= 0:
        raise ValueError("shell_type -1 is not supported.")
    return result


#
# gbasis wrappers
#


def _gob_cart_normalization(double alpha, long[::1] n not None):
    assert n.shape[0] == 3
    return gbasis.gob_cart_normalization(alpha, &n[0])


def _gob_pure_normalization(double alpha, long l):
    return gbasis.gob_pure_normalization(alpha, l)


cdef class GBasis:
    """
    This class describes basis sets applied to a certain molecular structure.

    The order of the pure shells is based on the order of real spherical.
    The functions are sorted from low to high magnetic quantum number,
    with cosine-like functions before the sine-like functions. The order
    of functions in a Cartesian shell is alphabetic. Some examples:

    shell_type = 0, S:
     0 -> 1
    shell_type = 1, P:
     0 -> x
     1 -> y
     2 -> z
    shell_type = 2, Cartesian D:
     0 -> xx
     1 -> xy
     2 -> xz
     3 -> yy
     4 -> yz
     5 -> zz
    shell_type = 3, Cartesian F:
     0 -> xxx
     1 -> xxy
     2 -> xxz
     3 -> xyy
     4 -> xyz
     5 -> xzz
     6 -> yyy
     7 -> yyz
     8 -> yzz
     9 -> zzz
    shell_type = -1, not allowed
    shell_type = -2, pure D:
     0 -> zz
     1 -> xz
     2 -> yz
     3 -> xx-yy
     4 -> xy
    shell_type = -3, pure F:
     0 -> zzz
     1 -> xzz
     2 -> yzz
     3 -> xxz-yyz
     4 -> xyz
     5 -> xxx-3xyy
     6 -> 3xxy-yyy

     Attributes
     ==========
     biblio : list
        A list of references to cite. It is populated depending on which methods were run during the
        calculation.
    """

    def __cinit__(self, centers: Iterable, shell_center: Iterable, shell_nprims: Iterable,
                  shell_types: Iterable, alphas: Iterable, con_coeffs: Iterable):
        """
        A C++ wrapper class for interfacing with all the gaussian integral and numerical integral
        code.

        Parameters
        ----------
        centers
            A numpy array with centers for the basis functions.
            shape = (ncenter, 3)
        shell_center
            An array with the center index for each shell.
            shape = (nshell,)
        shell_nprims
            The number of primitives in each shell.
            shape = (nshell,)
        shell_types
            An array with contraction types: 0 = S, 1 = P, 2 = Cartesian D,
            3 = Cartesian F, ..., -2 = pure D, -3 = pure F, ...
            shape = (nshell,)
        alphas
            The exponents of the primitives in one shell.
            shape = (sum(shell_nprims),)
        con_coeffs
            The contraction coefficients of the primitives for each
            contraction in a contiguous array. The coefficients are ordered
            according to the shells. Within each shell, the coefficients are
            grouped per exponent.
            shape = (sum(shell_nprims),)


        Copies are made of the arguments and stored internally, and are not meant
        to be modified once the GBasis object is created.

        Returns
        -------
        GBasis
            An instance of GBasis

        """
        # Make private copies of the input arrays.
        self._centers = np.array(centers, dtype=float)
        self._shell_center = np.array(shell_center, dtype=int)
        self._shell_nprims = np.array(shell_nprims, dtype=int)
        self._shell_types = np.array(shell_types, dtype=int)
        self._alphas = np.array(alphas, dtype=float)
        self._con_coeffs = np.array(con_coeffs, dtype=float)

        self._centers.flags.writeable = True
        # Set arrays unwritable because:
        #   (i) derived properties will be stored
        #   (ii) consistency tests below are only performed once (below)
        # In principle, one can make the arrays writable again, but then one
        # is deliberately taking risks. A strict read-only buffer would be
        # ideal but is not possible. One can always pass a pointer to a
        # C library and start messing things up.
        self._shell_center.flags.writeable = False
        self._shell_nprims.flags.writeable = False
        self._shell_types.flags.writeable = False
        self._alphas.flags.writeable = False
        self._con_coeffs.flags.writeable = False

        # Check array dimensions
        if self._centers.ndim != 2:
            raise TypeError('centers must be a 2D array')
        if self._shell_center.ndim != 1:
            raise TypeError('shell_center must be a 1D array')
        if self._shell_nprims.ndim != 1:
            raise TypeError('shell_nprims must be a 1D array')
        if self._shell_types.ndim != 1:
            raise TypeError('shell_types must be a 1D array')
        if self._alphas.ndim != 1:
            raise TypeError('alphas must be a 1D array')
        if self._con_coeffs.ndim != 1:
            raise TypeError('con_coeffs must be a 1D array')

        # Essential array checks
        if self._centers.shape[1] != 3:
            raise TypeError('centers must have three columns.')
        if self._shell_nprims.shape[0] != self._shell_center.shape[0]:
            raise TypeError('shell_nprims and shell_center must have the same length.')
        if self._shell_types.shape[0] != self._shell_center.shape[0]:
            raise TypeError('shell_types and shell_center must have the same length.')
        if self._alphas.shape[0] != self._con_coeffs.shape[0]:
            raise TypeError('alphas and con_coeffs must have the same length.')

        # Consistency checks
        if self._shell_center.min() < 0:
            raise ValueError('shell_center can not contain negative values.')
        if self._shell_center.max() >= self.centers.shape[0]:
            raise ValueError('shell_center can not contain values larger than the number of centers.')
        if self._shell_nprims.min() < 1:
            raise ValueError('shell_nprims elements must be strictly positive.')
        if (self._shell_types == -1).any():
            raise ValueError('The shell_type -1 is not supported.')
        cdef long nprim_total = self._shell_nprims.sum()
        if self._alphas.shape[0] != nprim_total:
            raise TypeError('The length of alphas must equal the total number of primitives.')


    def __init__(self, centers: Iterable, shell_center: Iterable, shell_nprims: Iterable,
                  shell_types: Iterable, alphas: Iterable, con_coeffs: Iterable):
        """
        A C++ wrapper class for interfacing with all the gaussian integral and numerical integral
        code.

        Parameters
        ----------
        centers
            A numpy array with centers for the basis functions.
            shape = (ncenter, 3)
        shell_center
            An array with the center index for each shell.
            shape = (nshell,)
        shell_nprims
            The number of primitives in each shell.
            shape = (nshell,)
        shell_types
            An array with contraction types: 0 = S, 1 = P, 2 = Cartesian D,
            3 = Cartesian F, ..., -2 = pure D, -3 = pure F, ...
            shape = (nshell,)
        alphas
            The exponents of the primitives in one shell.
            shape = (sum(shell_nprims),)
        con_coeffs
            The contraction coefficients of the primitives for each
            contraction in a contiguous array. The coefficients are ordered
            according to the shells. Within each shell, the coefficients are
            grouped per exponent.
            shape = (sum(shell_nprims),)


        Copies are made of the arguments and stored internally, and are not meant
        to be modified once the GOBasis object is created.

        Returns
        -------
        GBasis
            An instance of GOBasis

        """
        if self.__class__ == GBasis:
            raise NotImplementedError('GBasis is an abstract base class')
        self._log_init()
        self._biblio = []

    def __dealloc__(self):
        del self._baseptr

    @classmethod
    def concatenate(cls, *gbs: Type[GBasis]) -> Type[GBasis]:
        """Concatenate multiple basis objects into a new one.

        Parameters
        ----------
        gbs
            Each argument must be an instance of the same subclass of
            GBasis.

        Returns
        -------
        GBasis
            An instance of Gbasis (or its children) with its attributes concatenated together from
            the arguments.
        """

        # check if the classes match
        for gb in gbs:
            assert isinstance(gb, cls)

        # do the concatenation of each array properly
        centers = np.concatenate([gb.centers for gb in gbs])
        shell_center = []
        offset = 0
        for gb in gbs:
            shell_center.append(gb.shell_center + offset)
            offset += gb.ncenter
        shell_center = np.concatenate(shell_center)
        shell_nprims = np.concatenate([gb.shell_nprims for gb in gbs])
        shell_types = np.concatenate([gb.shell_types for gb in gbs])
        alphas = np.concatenate([gb.alphas for gb in gbs])
        con_coeffs = np.concatenate([gb.con_coeffs for gb in gbs])
        return cls(centers, shell_center, shell_nprims, shell_types, alphas, con_coeffs)

    @property
    def biblio(self):
        """References to cite. The list is generated depending on what users run."""
        return self._biblio

    # Array properties

    @property
    def centers(self):
        """The basis function centres"""
        return self._centers

    @property
    def shell_center(self):
        """The index in ``centers`` for each shell"""
        return self._shell_center.view()

    @property
    def shell_nprims(self):
        """The number of primitives in each shell."""
        return self._shell_nprims.view()

    @property
    def shell_types(self):
        """The type of each shell.
        0 = S, 1 = P, 2 = Cartesian D, 3 = Cartesian F, ..., -2 = pure D, -3 = pure F, ...
        """
        return self._shell_types.view()

    @property
    def alphas(self):
        """The exponents of the primitives in one shell."""
        return self._alphas.view()

    @property
    def con_coeffs(self):
        """The contraction coefficients of the primitives for each
        contraction in a contiguous array. The coefficients are ordered
        according to the shells. Within each shell, the coefficients are
        grouped per exponent.
        """
        return self._con_coeffs.view()

    # Array sizes

    @property
    def ncenter(self):
        """The number of centres"""
        return self.centers.shape[0]

    @property
    def nshell(self):
        """The number of shells"""
        return self.shell_center.shape[0]

    @property
    def nprim_total(self):
        """The total number of primitives"""
        return self.shell_nprims.sum()

    # Other properties

    @property
    def nbasis(self):
        """The number of basis functions"""
        return self._baseptr.get_nbasis()

    @property
    def shell_lookup(self):
        cdef np.npy_intp* shape = [self.nbasis]
        tmp = np.PyArray_SimpleNewFromData(1, shape,
                np.NPY_LONG, <void*> self._baseptr.get_shell_lookup())
        return tmp.copy()

    @property
    def basis_offsets(self):
        cdef np.npy_intp* shape = [self.nshell]
        tmp = np.PyArray_SimpleNewFromData(1, shape, np.NPY_LONG,
                    <void*> self._baseptr.get_basis_offsets())
        return tmp.copy()

    @property
    def nscales(self):
        """The number of normalization constants"""
        return self._baseptr.get_nscales()

    @property
    def max_shell_type(self):
        """The largest shell angular momentum in this basis object."""
        return self._baseptr.get_max_shell_type()

    def _log_init(self):
        """Write a summary of the basis to the screen logger"""
        # TODO: Re-enable log output
        # print('9: Initialized: %s' % self)
        # print(('9: Number of basis functions', self.nbasis))
        # print(('9: Number of normalization constants', self.nscales))
        # print(('9: Maximum shell type', self.max_shell_type))
        shell_type_names = {
            0: 'S', 1: 'P', 2: 'Dc', 3: 'Fc', 4:'Gc', 5: 'Hc', 6: 'Ic',
            -2: 'Dp', -3: 'Fp', -4:'Gp', -5: 'Hp', -6: 'Ip',
        }
        descs = ['']*self.ncenter
        for i in range(self.nshell):
            icenter = self.shell_center[i]
            s = descs[icenter]
            name = shell_type_names[self.shell_types[i]]
            s += ' %s%i' % (name, self.shell_nprims[i])
            descs[icenter] = s
        deflist = []
        # for i in range(self.ncenter):
        #     deflist.append(('9: Center % 5i' % i, descs[i]))
        # print(deflist)
        # print()

    def get_scales(self):
        """Return a copy of the normalization constants."""
        cdef np.npy_intp shape[1]
        shape[0] = self.nscales
        tmp = np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, <void*> self._baseptr.get_scales(0))
        return tmp.copy()

    # low-level compute routines, for debugging only
    def _compute_grid_point1(self, double[::1] output not None,
                            double[::1] point not None,
                            _GB1DMGridFn grid_fn not None):
        assert output.shape[0] == self.nbasis
        assert point.shape[0] == 3
        self._baseptr.compute_grid_point1(&output[0], &point[0], grid_fn._baseptr)

    def get_subset(self, ishells: List) -> Type[GBasis]:
        """Construct a sub basis set for a selection of shells

        Parameters
        ----------
        ishells
            A list of indexes of shells to be retained in the sub basis set

        Returns
        -------
        GBasis
            An instance of the same class as self containing only
            the basis functions that correspond to the select shells in
            the ``ishells`` list.
        """
        # find the centers corresponding to ishells
        icenters = set([])
        for ishell in ishells:
            if ishell < 0:
                raise ValueError('ishell out of range: %i < 0' % ishell)
            if ishell >= self.nshell:
                raise ValueError('ishell out of range: %i >= %s' % (ishell, self.nshell))
            icenters.add(self.shell_center[ishell])
        icenters = sorted(icenters) # fix the order
        new_centers = self.centers[icenters]
        # construct the new shell_center, shell_nprims, shell_types
        new_shell_center = np.zeros(len(ishells), int)
        new_nprims = np.zeros(len(ishells), int)
        new_shell_types = np.zeros(len(ishells), int)
        for new_ishell, ishell in enumerate(ishells):
            new_shell_center[new_ishell] = icenters.index(self.shell_center[ishell])
            new_nprims[new_ishell] = self.shell_nprims[ishell]
            new_shell_types[new_ishell] = self.shell_types[ishell]
        # construct the new alphas and con_coeffs
        new_nprim_total = new_nprims.sum()
        new_alphas = np.zeros(new_nprim_total, float)
        new_con_coeffs = np.zeros(new_nprim_total, float)
        new_iprim = 0
        for new_ishell, ishell in enumerate(ishells):
            nprim = new_nprims[new_ishell]
            old_iprim = self.shell_nprims[:ishell].sum()
            new_alphas[new_iprim:new_iprim+nprim] = self.alphas[old_iprim:old_iprim+nprim]
            new_con_coeffs[new_iprim:new_iprim+nprim] = self.con_coeffs[old_iprim:old_iprim+nprim]
            new_iprim += nprim
        # make a mapping between the indices of old and new basis functions
        ibasis_list = []
        for new_ishell, ishell in enumerate(ishells):
            ibasis_old = sum(common.get_shell_nbasis(self.shell_types[i]) for i in range(ishell))
            nbasis = common.get_shell_nbasis(self.shell_types[ishell])
            ibasis_list.extend(range(ibasis_old, ibasis_old+nbasis))
        ibasis_list = np.array(ibasis_list)
        # create the basis set object
        basis = self.__class__(new_centers, new_shell_center, new_nprims,
                               new_shell_types, new_alphas, new_con_coeffs)
        # return stuff
        return basis, ibasis_list

    def get_basis_atoms(self, coordinates: np.ndarray) -> List[Tuple[GBasis, List]]:
        """Return a list of atomic basis sets for a set of coordinates

        Parameters
        ----------
        coordinates
            An (N, 3) array with atomic coordinates, used to find the
            centers associated with atoms. An exact match of the Cartesian
            coordinates is required to properly select a shell.

        Returns
        -------
        List
            A list with one tuple for every atom: (gbasis,
            ibasis_list), where gbasis is a basis set object for the atom and
            ibasis_list is a list of basis set indexes that can be used to
            substitute results from the atomic basis set back into the molecular
            basis set.

            For example, when a density matrix for the atom is
            obtained and it needs to be plugged back into the molecular density
            matrix, one can do the following::

            >>>> mol_dm[ibasis_list, ibasis_list.reshape(-1,1)] = atom_dm

        """
        result = []
        for c in coordinates:
            # find the corresponding center(s).
            icenters = []
            for icenter in range(self.ncenter):
                # require an exact match of the coordinates
                if (self.centers[icenter] == c).all():
                    icenters.append(icenter)
            icenters = set(icenters)
            # find the shells on these centers
            ishells = []
            for ishell in range(self.nshell):
                if self.shell_center[ishell] in icenters:
                    ishells.append(ishell)
            # construct a sub basis
            sub_basis, ibasis_list = self.get_subset(ishells)
            result.append((sub_basis, ibasis_list))
        return result


cdef class GOBasis(GBasis):
    def __cinit__(self, centers, shell_center, shell_nprims, shell_types, alphas, con_coeffs):
        self._this = new gbasis.GOBasis(
            <double*>self._centers.data, <long*>self._shell_center.data,
            <long*>self._shell_nprims.data, <long*>self._shell_types.data,
            <double*>self._alphas.data, <double*>self._con_coeffs.data,
            self._centers.shape[0], self._shell_types.shape[0], self._alphas.shape[0])
        self._baseptr = <gbasis.GBasis*> self._this

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


    def _compute_cholesky(self, _GB4Integral gb4int, double threshold=1e-8) -> np.ndarray:
        """Apply the Cholesky code to a given type of four-center integrals.

        Parameters
        ----------
        gb4int
            The object that can carry out four-center integrals.
        threshold
            The cutoff for the Cholesky decomposition.

        Returns
        -------
        np.ndarray, shape(nvec, nbasis, nbasis), dtype=float
            The Cholesky-decomposed four-center integrals
        """
        cdef gbw.GB4IntegralWrapper* gb4w = NULL
        cdef vector[double]* vectors = NULL
        cdef np.npy_intp dims[3]
        cdef np.ndarray result

        try:
            gb4w = new gbw.GB4IntegralWrapper(self._this,
                                              gb4int._baseptr)
            vectors = new vector[double]()
            nvec = cholesky.cholesky(gb4w, vectors, threshold)
            dims[0] = <np.npy_intp> nvec
            dims[1] = <np.npy_intp> self.nbasis
            dims[2] = <np.npy_intp> self.nbasis
            result = np.PyArray_SimpleNewFromData(3, dims, np.NPY_DOUBLE, &(deref(vectors)[0]))
            PyArray_ENABLEFLAGS(result, np.NPY_OWNDATA)
        finally:
            if gb4w is not NULL:
                del gb4w

        return result

    def compute_electron_repulsion_cholesky(self, double threshold=1e-8) -> np.ndarray:
        r"""Compute Cholesky decomposition of electron repulsion four-center integrals.

        .. math::
            L_{ik}^{v}\cdot L_{jl}^{v}=\mel{\chi_{i}\chi_{j}}{\frac{1}{\abs{r_{1}-r_{2}}}}{\chi_{k}\chi_{l}}

        Parameters
        ----------
        threshold
            The cutoff for the Cholesky decomposition.

        Returns
        -------
        np.ndarray, shape(nvec, nbasis, nbasis), dtype=float
            The Cholesky-decomposed four-center integrals

        """
        return self._compute_cholesky(_GB4ElectronRepulsionIntegralLibInt(self.max_shell_type))

    def compute_erf_repulsion_cholesky(self, double mu=0.0, double threshold=1e-8) -> np.ndarray:
        r"""Compute Cholesky decomposition of Erf repulsion four-center integrals.

        .. math::
            L_{ik}^{v}\cdot L_{jl}^{v}=\mel{\chi_{i}\chi_{j}}{\frac{\mathrm{erf}(\mu r)}{r}}{\chi_{k}\chi_{l}}

        Parameters
        ----------
        mu : float
            Parameter for the erf(mu r)/r potential. Default is zero.
        threshold
            The cutoff for the Cholesky decomposition.

        Returns
        -------
        np.ndarray, shape(nvec, nbasis, nbasis), dtype=float
            The Cholesky-decomposed four-center integrals

        """
        return self._compute_cholesky(_GB4ErfIntegralLibInt(self.max_shell_type, mu))

    def compute_gauss_repulsion_cholesky(self, double c=1.0, double alpha=1.0,
                                         double threshold=1e-8) -> np.ndarray:
        r"""Compute Cholesky decomposition of Gauss repulsion four-center integrals.

        .. math::
            L_{ik}^{v}\cdot L_{jl}^{v}=\mel{\chi_{i}\chi_{j}}{c\exp(-\alpha r^{2})}{\chi_{k}\chi_{l}}

        Parameters
        ----------
        c : float
            Coefficient of the gaussian.
        alpha : float
            Exponential parameter of the gaussian.
        threshold
            The cutoff for the Cholesky decomposition.

        Returns
        -------
        np.ndarray, shape(nvec, nbasis, nbasis), dtype=float
            The Cholesky-decomposed four-center integrals


        """
        return self._compute_cholesky(_GB4GaussIntegralLibInt(self.max_shell_type, c, alpha))

    def compute_ralpha_repulsion_cholesky(self, double alpha=-1.0,
                                          double threshold=1e-8) -> np.ndarray:
        r"""Compute Cholesky decomposition of ralpha repulsion four-center integrals.

        .. math::
            L_{ik}^{v}\cdot L_{jl}^{v}=\mel{\chi_{i}\chi_{j}}{r^{\alpha}}{\chi_{k}\chi_{l}}

        with :math:`\alpha > -3`.
        Parameters
        ----------
        alpha : float
            The power of r in the interaction potential.
        threshold
            The cutoff for the Cholesky decomposition.

        Returns
        -------
        np.ndarray, shape(nvec, nbasis, nbasis), dtype=float
            The Cholesky-decomposed four-center integrals


        """
        return self._compute_cholesky(_GB4RAlphaIntegralLibInt(self.max_shell_type, alpha))

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
        pot_stride = (pots.strides[0]/8)
        if pots.shape[1] > 1:
            pot_stride *= (pots.strides[1]/8)
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
# ints wrappers (for testing and use in Cholesky iterators only)
#

ints.libint2_static_init()
def libint2_static_cleanup():
    ints.libint2_static_cleanup()
atexit.register(libint2_static_cleanup)


cdef class _GB4Integral:
    """Wrapper for ints.GB4Integral. For testing only."""
    cdef ints.GB4Integral* _baseptr

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

    def reset(self, long shell_type0, long shell_type1, long shell_type2, long shell_type3,
              double[::1] r0 not None, double[::1] r1 not None,
              double[::1] r2 not None, double[::1] r3 not None):
        assert r0.shape[0] == 3
        assert r1.shape[0] == 3
        assert r2.shape[0] == 3
        assert r3.shape[0] == 3
        self._baseptr.reset(shell_type0, shell_type1, shell_type2, shell_type3,
                         &r0[0], &r1[0], &r2[0], &r3[0])

    def add(self, double coeff, double alpha0, double alpha1, double alpha2, double alpha3,
            double[::1] scales0 not None, double[::1] scales1 not None,
            double[::1] scales2 not None, double[::1] scales3 not None):
        assert scales0.shape[0] == _get_shell_nbasis(abs(self._baseptr.get_shell_type0()))
        assert scales1.shape[0] == _get_shell_nbasis(abs(self._baseptr.get_shell_type1()))
        assert scales2.shape[0] == _get_shell_nbasis(abs(self._baseptr.get_shell_type2()))
        assert scales3.shape[0] == _get_shell_nbasis(abs(self._baseptr.get_shell_type3()))
        self._baseptr.add(coeff, alpha0, alpha1, alpha2, alpha3,
                       &scales0[0], &scales1[0],
                       &scales2[0], &scales3[0])

    def cart_to_pure(self):
        self._baseptr.cart_to_pure()

    def get_work(self, shape0, shape1, shape2, shape3):
        """This returns a **copy** of the c++ work array.

           Returning a numpy array with a buffer created in c++ is dangerous.
           If the c++ array becomes deallocated, the numpy array may still
           point to the deallocated memory. For that reason, a copy is returned.
           Speed is not an issue as this class is only used for testing.
        """
        cdef np.npy_intp shape[4]
        assert shape0 > 0
        assert shape1 > 0
        assert shape2 > 0
        assert shape3 > 0
        assert shape0 <= self.max_nbasis
        assert shape1 <= self.max_nbasis
        assert shape2 <= self.max_nbasis
        assert shape3 <= self.max_nbasis
        shape[0] = shape0
        shape[1] = shape1
        shape[2] = shape2
        shape[3] = shape3
        tmp = np.PyArray_SimpleNewFromData(4, shape, np.NPY_DOUBLE, <void*> self._baseptr.get_work())
        return tmp.copy()


cdef class _GB4ElectronRepulsionIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4ElectronRepulsionIntegralLibInt, for testing only"""
    cdef ints.GB4ElectronRepulsionIntegralLibInt* _this

    def __cinit__(self, long max_nbasis):
        self._this = new ints.GB4ElectronRepulsionIntegralLibInt(max_nbasis)
        self._baseptr = <ints.GB4Integral*> self._this


cdef class _GB4ErfIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4ElectronRepulsionIntegralLibInt, for testing only"""
    cdef ints.GB4ErfIntegralLibInt* _this

    def __cinit__(self, long max_nbasis, double mu):
        self._this = new ints.GB4ErfIntegralLibInt(max_nbasis, mu)
        self._baseptr = <ints.GB4Integral*> self._this

    @property
    def mu(self):
        return self._this.get_mu()


cdef class _GB4GaussIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4GaussIntegralLibInt, for testing only"""
    cdef ints.GB4GaussIntegralLibInt* _this

    def __cinit__(self, long max_nbasis, double c, double alpha):

        self._this = new ints.GB4GaussIntegralLibInt(max_nbasis, c, alpha)
        self._baseptr = <ints.GB4Integral*> self._this

    @property
    def c(self):
        return self._this.get_c()

    @property
    def alpha(self):
        return self._this.get_alpha()


cdef class _GB4RAlphaIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4RAlphaIntegralLibInt, for testing only"""
    cdef ints.GB4RAlphaIntegralLibInt* _this

    def __cinit__(self, long max_nbasis, double alpha):
        self._this = new ints.GB4RAlphaIntegralLibInt(max_nbasis, alpha)
        self._baseptr = <ints.GB4Integral*> self._this

    @property
    def alpha(self):
        return self._this.get_alpha()

cdef class _GB4DeltaIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4DeltaIntegralLibInt, for testing only"""
    cdef ints.GB4DeltaIntegralLibInt* _this

    def __cinit__(self, long max_nbasis):
        self._this = new ints.GB4DeltaIntegralLibInt(max_nbasis)
        self._baseptr = <ints.GB4Integral*> self._this


cdef class _GB4IntraDensIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4IntraDensIntegralLibInt, for testing only"""
    cdef ints.GB4IntraDensIntegralLibInt* _this

    def __cinit__(self, long max_nbasis, double[:, ::1] point not None):
        self._this = new ints.GB4IntraDensIntegralLibInt(max_nbasis, &point[0, 0])
        self._baseptr = <ints.GB4Integral*> self._this



#
# fns wrappers (for testing and use in this module)
#


cdef class _GB1DMGridFn:
    """Wrapper for fns.GB1DMGridFn, for testing only"""
    cdef fns.GB1DMGridFn* _baseptr

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
    cdef fns.GB1DMGridDensityFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridDensityFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridGradientFn(_GB1DMGridFn):
    cdef fns.GB1DMGridGradientFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridGradientFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridGGAFn(_GB1DMGridFn):
    cdef fns.GB1DMGridGGAFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridGGAFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridKineticFn(_GB1DMGridFn):
    cdef fns.GB1DMGridKineticFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridKineticFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridHessianFn(_GB1DMGridFn):
    cdef fns.GB1DMGridHessianFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridHessianFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridMGGAFn(_GB1DMGridFn):
    cdef fns.GB1DMGridMGGAFn* _this
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridMGGAFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this








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
    nucpot.compute_grid_nucpot(
        &coordinates[0,0], &charges[0], ncharge,
        &points[0,0], &output[0], npoint)
