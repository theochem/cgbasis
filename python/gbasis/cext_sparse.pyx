#cython: language_level=3

import atexit

cimport numpy as np
import numpy as np

np.import_array()

from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

from gbasis.pxds.c_numpy_wrapper cimport PyArray_ENABLEFLAGS

from gbasis.cext_common cimport GBasis, _check_shape, _prepare_array, _get_shell_nbasis
from gbasis.pxds cimport c_gbasis
from gbasis.pxds.sparse cimport c_gbw
from gbasis.pxds.sparse cimport c_cholesky
from gbasis.pxds.twos cimport c_ints4

cdef class GOBasisSparse(GBasis):
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
        cdef c_gbw.GB4IntegralWrapper* gb4w = NULL
        cdef vector[double]* vectors = NULL
        cdef np.npy_intp dims[3]
        cdef np.ndarray result

        try:
            gb4w = new c_gbw.GB4IntegralWrapper(self._this,
                                              gb4int._baseptr)
            vectors = new vector[double]()
            nvec = c_cholesky.cholesky(gb4w, vectors, threshold)
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

#
# ints wrappers (for testing and use in Cholesky iterators only)
#

c_ints4.libint2_static_init()
def libint2_static_cleanup():
    c_ints4.libint2_static_cleanup()
atexit.register(libint2_static_cleanup)


cdef class _GB4Integral:
    """Wrapper for c_ints4.GB4Integral. For testing only."""
    cdef c_ints4.GB4Integral* _baseptr

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
    """Wrapper for c_ints4.GB4ElectronRepulsionIntegralLibInt, for testing only"""
    cdef c_ints4.GB4ElectronRepulsionIntegralLibInt* _this

    def __cinit__(self, long max_nbasis):
        self._this = new c_ints4.GB4ElectronRepulsionIntegralLibInt(max_nbasis)
        self._baseptr = <c_ints4.GB4Integral*> self._this


cdef class _GB4ErfIntegralLibInt(_GB4Integral):
    """Wrapper for c_ints4.GB4ElectronRepulsionIntegralLibInt, for testing only"""
    cdef c_ints4.GB4ErfIntegralLibInt* _this

    def __cinit__(self, long max_nbasis, double mu):
        self._this = new c_ints4.GB4ErfIntegralLibInt(max_nbasis, mu)
        self._baseptr = <c_ints4.GB4Integral*> self._this

    @property
    def mu(self):
        return self._this.get_mu()


cdef class _GB4GaussIntegralLibInt(_GB4Integral):
    """Wrapper for c_ints4.GB4GaussIntegralLibInt, for testing only"""
    cdef c_ints4.GB4GaussIntegralLibInt* _this

    def __cinit__(self, long max_nbasis, double c, double alpha):

        self._this = new c_ints4.GB4GaussIntegralLibInt(max_nbasis, c, alpha)
        self._baseptr = <c_ints4.GB4Integral*> self._this

    @property
    def c(self):
        return self._this.get_c()

    @property
    def alpha(self):
        return self._this.get_alpha()


cdef class _GB4RAlphaIntegralLibInt(_GB4Integral):
    """Wrapper for c_ints4.GB4RAlphaIntegralLibInt, for testing only"""
    cdef c_ints4.GB4RAlphaIntegralLibInt* _this

    def __cinit__(self, long max_nbasis, double alpha):
        self._this = new c_ints4.GB4RAlphaIntegralLibInt(max_nbasis, alpha)
        self._baseptr = <c_ints4.GB4Integral*> self._this

    @property
    def alpha(self):
        return self._this.get_alpha()

cdef class _GB4DeltaIntegralLibInt(_GB4Integral):
    """Wrapper for c_ints4.GB4DeltaIntegralLibInt, for testing only"""
    cdef c_ints4.GB4DeltaIntegralLibInt* _this

    def __cinit__(self, long max_nbasis):
        self._this = new c_ints4.GB4DeltaIntegralLibInt(max_nbasis)
        self._baseptr = <c_ints4.GB4Integral*> self._this


cdef class _GB4IntraDensIntegralLibInt(_GB4Integral):
    """Wrapper for c_ints4.GB4IntraDensIntegralLibInt, for testing only"""
    cdef c_ints4.GB4IntraDensIntegralLibInt* _this

    def __cinit__(self, long max_nbasis, double[:, ::1] point not None):
        self._this = new c_ints4.GB4IntraDensIntegralLibInt(max_nbasis, &point[0, 0])
        self._baseptr = <c_ints4.GB4Integral*> self._this
