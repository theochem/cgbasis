# cython: embedsignature=True, language_level=3

from typing import Iterable, Tuple, Type, List, Union

cimport numpy as np
import numpy as np

np.import_array()

cimport c_common
cimport c_gbasis


cdef _check_shape(array, shape, name):
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


cdef _prepare_array(array, shape, name):
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
    result = c_common.get_shell_nbasis(shell_type)
    if result <= 0:
        raise ValueError("shell_type -1 is not supported.")
    return result


#
# gbasis wrappers
#


def _gob_cart_normalization(double alpha, long[::1] n not None):
    assert n.shape[0] == 3
    return c_gbasis.gob_cart_normalization(alpha, &n[0])


def _gob_pure_normalization(double alpha, long l):
    return c_gbasis.gob_pure_normalization(alpha, l)


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

    # cdef c_gbasis.GBasis* _baseptr
    # # Keep reference to arrays to make sure they will not be deallocated.
    # cdef np.ndarray _centers
    # cdef np.ndarray _shell_map
    # cdef np.ndarray _nprims
    # cdef np.ndarray _shell_types
    # cdef np.ndarray _alphas
    # cdef np.ndarray _con_coeffs

    def __cinit__(self, centers: Iterable, shell_map: Iterable, nprims: Iterable,
                  shell_types: Iterable, alphas: Iterable, con_coeffs: Iterable):
        """
        A C++ wrapper class for interfacing with all the gaussian integral and numerical integral
        code.

        Parameters
        ----------
        centers
            A numpy array with centers for the basis functions.
            shape = (ncenter, 3)
        shell_map
            An array with the center index for each shell.
            shape = (nshell,)
        nprims
            The number of primitives in each shell.
            shape = (nshell,)
        shell_types
            An array with contraction types: 0 = S, 1 = P, 2 = Cartesian D,
            3 = Cartesian F, ..., -2 = pure D, -3 = pure F, ...
            shape = (nshell,)
        alphas
            The exponents of the primitives in one shell.
            shape = (sum(nprims),)
        con_coeffs
            The contraction coefficients of the primitives for each
            contraction in a contiguous array. The coefficients are ordered
            according to the shells. Within each shell, the coefficients are
            grouped per exponent.
            shape = (sum(nprims),)


        Copies are made of the arguments and stored internally, and are not meant
        to be modified once the GBasis object is created.

        Returns
        -------
        GBasis
            An instance of GBasis

        """
        # Make private copies of the input arrays.
        self._centers = np.array(centers, dtype=float)
        self._shell_map = np.array(shell_map, dtype=int)
        self._nprims = np.array(nprims, dtype=int)
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
        self._shell_map.flags.writeable = False
        self._nprims.flags.writeable = False
        self._shell_types.flags.writeable = False
        self._alphas.flags.writeable = False
        self._con_coeffs.flags.writeable = False

        # Check array dimensions
        if self._centers.ndim != 2:
            raise TypeError('centers must be a 2D array')
        if self._shell_map.ndim != 1:
            raise TypeError('shell_map must be a 1D array')
        if self._nprims.ndim != 1:
            raise TypeError('nprims must be a 1D array')
        if self._shell_types.ndim != 1:
            raise TypeError('shell_types must be a 1D array')
        if self._alphas.ndim != 1:
            raise TypeError('alphas must be a 1D array')
        if self._con_coeffs.ndim != 1:
            raise TypeError('con_coeffs must be a 1D array')

        # Essential array checks
        if self._centers.shape[1] != 3:
            raise TypeError('centers must have three columns.')
        if self._nprims.shape[0] != self._shell_map.shape[0]:
            raise TypeError('nprims and shell_map must have the same length.')
        if self._shell_types.shape[0] != self._shell_map.shape[0]:
            raise TypeError('shell_types and shell_map must have the same length.')
        if self._alphas.shape[0] != self._con_coeffs.shape[0]:
            raise TypeError('alphas and con_coeffs must have the same length.')

        # Consistency checks
        if self._shell_map.min() < 0:
            raise ValueError('shell_map can not contain negative values.')
        if self._shell_map.max() >= self.centers.shape[0]:
            raise ValueError('shell_map can not contain values larger than the number of centers.')
        if self._nprims.min() < 1:
            raise ValueError('nprims elements must be strictly positive.')
        if (self._shell_types == -1).any():
            raise ValueError('The shell_type -1 is not supported.')
        cdef long nprim_total = self._nprims.sum()
        if self._alphas.shape[0] != nprim_total:
            raise TypeError('The length of alphas must equal the total number of primitives.')




    def __init__(self, centers: Iterable, shell_map: Iterable, nprims: Iterable,
                  shell_types: Iterable, alphas: Iterable, con_coeffs: Iterable):
        """
        A C++ wrapper class for interfacing with all the gaussian integral and numerical integral
        code.

        Parameters
        ----------
        centers
            A numpy array with centers for the basis functions.
            shape = (ncenter, 3)
        shell_map
            An array with the center index for each shell.
            shape = (nshell,)
        nprims
            The number of primitives in each shell.
            shape = (nshell,)
        shell_types
            An array with contraction types: 0 = S, 1 = P, 2 = Cartesian D,
            3 = Cartesian F, ..., -2 = pure D, -3 = pure F, ...
            shape = (nshell,)
        alphas
            The exponents of the primitives in one shell.
            shape = (sum(nprims),)
        con_coeffs
            The contraction coefficients of the primitives for each
            contraction in a contiguous array. The coefficients are ordered
            according to the shells. Within each shell, the coefficients are
            grouped per exponent.
            shape = (sum(nprims),)


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
        shell_map = []
        offset = 0
        for gb in gbs:
            shell_map.append(gb.shell_map + offset)
            offset += gb.ncenter
        shell_map = np.concatenate(shell_map)
        nprims = np.concatenate([gb.nprims for gb in gbs])
        shell_types = np.concatenate([gb.shell_types for gb in gbs])
        alphas = np.concatenate([gb.alphas for gb in gbs])
        con_coeffs = np.concatenate([gb.con_coeffs for gb in gbs])
        return cls(centers, shell_map, nprims, shell_types, alphas, con_coeffs)

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
    def shell_map(self):
        """The index in ``centers`` for each shell"""
        return self._shell_map.view()

    @property
    def nprims(self):
        """The number of primitives in each shell."""
        return self._nprims.view()

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
        return self.shell_map.shape[0]

    @property
    def nprim_total(self):
        """The total number of primitives"""
        return self.nprims.sum()

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
            icenter = self.shell_map[i]
            s = descs[icenter]
            name = shell_type_names[self.shell_types[i]]
            s += ' %s%i' % (name, self.nprims[i])
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
            icenters.add(self.shell_map[ishell])
        icenters = sorted(icenters) # fix the order
        new_centers = self.centers[icenters]
        # construct the new shell_map, nprims, shell_types
        new_shell_map = np.zeros(len(ishells), int)
        new_nprims = np.zeros(len(ishells), int)
        new_shell_types = np.zeros(len(ishells), int)
        for new_ishell, ishell in enumerate(ishells):
            new_shell_map[new_ishell] = icenters.index(self.shell_map[ishell])
            new_nprims[new_ishell] = self.nprims[ishell]
            new_shell_types[new_ishell] = self.shell_types[ishell]
        # construct the new alphas and con_coeffs
        new_nprim_total = new_nprims.sum()
        new_alphas = np.zeros(new_nprim_total, float)
        new_con_coeffs = np.zeros(new_nprim_total, float)
        new_iprim = 0
        for new_ishell, ishell in enumerate(ishells):
            nprim = new_nprims[new_ishell]
            old_iprim = self.nprims[:ishell].sum()
            new_alphas[new_iprim:new_iprim+nprim] = self.alphas[old_iprim:old_iprim+nprim]
            new_con_coeffs[new_iprim:new_iprim+nprim] = self.con_coeffs[old_iprim:old_iprim+nprim]
            new_iprim += nprim
        # make a mapping between the indices of old and new basis functions
        ibasis_list = []
        for new_ishell, ishell in enumerate(ishells):
            ibasis_old = sum(c_common.get_shell_nbasis(self.shell_types[i]) for i in range(ishell))
            nbasis = c_common.get_shell_nbasis(self.shell_types[ishell])
            ibasis_list.extend(range(ibasis_old, ibasis_old+nbasis))
        ibasis_list = np.array(ibasis_list)
        # create the basis set object
        basis = self.__class__(new_centers, new_shell_map, new_nprims,
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
                if self.shell_map[ishell] in icenters:
                    ishells.append(ishell)
            # construct a sub basis
            sub_basis, ibasis_list = self.get_subset(ishells)
            result.append((sub_basis, ibasis_list))
        return result
