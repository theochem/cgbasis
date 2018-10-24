cimport numpy as np
import numpy as np

np.import_array() # needed when cimporting np

cimport fns
from gbasis.cext cimport _get_shell_nbasis

__all__ = [
    # fns
    '_GB1DMGridDensityFn', '_GB1DMGridGradientFn', '_GB1DMGridGGAFn',
    '_GB1DMGridKineticFn', '_GB1DMGridHessianFn', '_GB1DMGridMGGAFn',
]

#
# fns wrappers (for testing and use in cext.pyx)
#


cdef class _GB1DMGridFn:
    """Wrapper for fns.GB1DMGridFn, for testing only"""

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
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridDensityFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridGradientFn(_GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridGradientFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridGGAFn(_GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridGGAFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridKineticFn(_GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridKineticFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridHessianFn(_GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridHessianFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this


cdef class _GB1DMGridMGGAFn(_GB1DMGridFn):
    def __cinit__(self, long max_nbasis):
        self._this = new fns.GB1DMGridMGGAFn(max_nbasis)
        self._baseptr = <fns.GB1DMGridFn*> self._this