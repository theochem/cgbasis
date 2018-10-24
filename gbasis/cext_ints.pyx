cimport numpy as np
import numpy as np

np.import_array() # needed when cimporting np

cimport ints

from gbasis.cext cimport _get_shell_nbasis

import atexit

__all__ = [
    # ints
    '_GB4Integral',
    '_GB4ElectronRepulsionIntegralLibInt',
    '_GB4ErfIntegralLibInt', '_GB4GaussIntegralLibInt',
    '_GB4RAlphaIntegralLibInt',
    '_GB4DeltaIntegralLibInt',
    '_GB4IntraDensIntegralLibInt',
]

#
# ints wrappers (for testing and use in Cholesky iterators only)
#

ints.libint2_static_init()
def libint2_static_cleanup():
    ints.libint2_static_cleanup()
atexit.register(libint2_static_cleanup)


cdef class _GB4Integral:
    """Wrapper for ints.GB4Integral. For testing only."""

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

    def __cinit__(self, long max_nbasis):
        self._this = new ints.GB4ElectronRepulsionIntegralLibInt(max_nbasis)
        self._baseptr = <ints.GB4Integral*> self._this


cdef class _GB4ErfIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4ElectronRepulsionIntegralLibInt, for testing only"""

    def __cinit__(self, long max_nbasis, double mu):
        self._this = new ints.GB4ErfIntegralLibInt(max_nbasis, mu)
        self._baseptr = <ints.GB4Integral*> self._this

    @property
    def mu(self):
        return self._this.get_mu()


cdef class _GB4GaussIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4GaussIntegralLibInt, for testing only"""

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

    def __cinit__(self, long max_nbasis, double alpha):
        self._this = new ints.GB4RAlphaIntegralLibInt(max_nbasis, alpha)
        self._baseptr = <ints.GB4Integral*> self._this

    @property
    def alpha(self):
        return self._this.get_alpha()

cdef class _GB4DeltaIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4DeltaIntegralLibInt, for testing only"""

    def __cinit__(self, long max_nbasis):
        self._this = new ints.GB4DeltaIntegralLibInt(max_nbasis)
        self._baseptr = <ints.GB4Integral*> self._this


cdef class _GB4IntraDensIntegralLibInt(_GB4Integral):
    """Wrapper for ints.GB4IntraDensIntegralLibInt, for testing only"""

    def __cinit__(self, long max_nbasis, double[:, ::1] point not None):
        self._this = new ints.GB4IntraDensIntegralLibInt(max_nbasis, &point[0, 0])
        self._baseptr = <ints.GB4Integral*> self._this