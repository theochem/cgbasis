cimport numpy as np
np.import_array()

import numpy as np # Must include this line after cimport numpy

cimport libc.string

from gbasis cimport iter_pow
from gbasis cimport iter_gb
from gbasis cimport ints
from gbasis cimport gbw
from gbasis cimport common
from gbasis cimport boys
from gbasis cimport cartpure
from gbasis.cext cimport GBasis, GOBasis, _get_shell_nbasis

__all__ = [
    # boys
    '_boys_function', '_boys_function_array',
    # cartpure
    '_cart_to_pure_low',
    # common.cpp
    '_fac', '_fac2', '_binom',
    '_get_max_shell_type',
    '_gpt_coeff', '_gb_overlap_int1d', '_nuclear_attraction_helper',
    '_cit', '_jfac', '_dtaylor',
    # gbw
    '_get_2index_slice', '_compute_diagonal', '_select_2index',
    # ints
    '_GB2OverlapIntegral', '_GB2KineticIntegral',
    '_GB2ErfAttractionIntegral',
    '_GB2GaussAttractionIntegral',
    '_GB2NuclearAttractionIntegral',
    # iter_gb
    '_IterGB1', '_IterGB2', '_IterGB4',
    # iter_pow
    '_iter_pow1_inc', '_IterPow1', '_IterPow2',
]

#
# boys wrappers (for testing only)
#


def _boys_function(long m, double t):
    return boys.boys_function(m, t)


def _boys_function_array(long mmax, double t):
    cdef double[::1] output = np.zeros(mmax+1)
    boys.boys_function_array(mmax, t, &output[0])
    return output


#
# cartpure wrappers (for testing only)
#


def _cart_to_pure_low(double[::1] work_cart not None,
                     double[::1] work_pure not None,
                     long shell_type, long nant, long npost):
    cartpure.cart_to_pure_low(
        &work_cart[0], &work_pure[0], shell_type, nant,
        npost
    )

#
# common.cpp wrappers
#

def _fac(long n):
    return common.fac(n)


def _fac2(long n):
    return common.fac2(n)


def _binom(long n, long m):
    return common.binom(n, m)

def _get_max_shell_type():
    return common.get_max_shell_type()


def _gpt_coeff(long k, long n0, long n1, double pa, double pb):
    return common.gpt_coeff(k, n0, n1, pa, pb)


def _gb_overlap_int1d(long n0, long n1, double pa, double pb, double inv_gamma):
    return common.gb_overlap_int1d(n0, n1, pa, pb, inv_gamma)


def _nuclear_attraction_helper(double[::1] work_g not None,
                                long n0, long n1, double pa, double pb, double cp,
                                double gamma_inv):
    assert work_g.shape[0] == n0+n1+1
    common.nuclear_attraction_helper(&work_g[0], n0, n1, pa, pb, cp, gamma_inv)

def _cit(int i, double t, int m):
        return common.cit(i, t, m)

def _jfac(int j, int n):
    return common.jfac(j, n)

def _dtaylor(int n, double alpha, double t, double tfactor):
    return common.dtaylor(n, alpha, t, tfactor)


#
# gbw wrappers
#

def _select_2index(GOBasis gobasis, long index0, long index2):
    """Select 2index 2e ints to compute. For cholesky testing only"""
    assert 0 <= index0 < gobasis.nbasis
    assert 0 <= index2 < gobasis.nbasis

    cdef ints.GB4ElectronRepulsionIntegralLibInt* gb4int = NULL
    cdef gbw.GB4IntegralWrapper* gb4w = NULL

    cdef long pbegin0
    cdef long pend0
    cdef long pbegin2
    cdef long pend2

    try:
        gb4int = new ints.GB4ElectronRepulsionIntegralLibInt(
                            gobasis.max_shell_type)
        gb4w = new gbw.GB4IntegralWrapper(gobasis._this,
                            <ints.GB4Integral*> gb4int)
        gb4w.select_2index(index0, index2, &pbegin0, &pend0, &pbegin2, &pend2)
    finally:
        if gb4int is not NULL:
            del gb4int
        if gb4w is not NULL:
            del gb4w
    return pbegin0, pend0, pbegin2, pend2


def _compute_diagonal(GOBasis gobasis, double[:, ::1] diagonal not None):
    """Get the diagonal 2e ints. For cholesky testing only"""
    cdef ints.GB4ElectronRepulsionIntegralLibInt* gb4int = NULL
    cdef gbw.GB4IntegralWrapper* gb4w = NULL
    cdef double[:, ::1] output
    output = diagonal

    try:
        gb4int = new ints.GB4ElectronRepulsionIntegralLibInt(
                            gobasis.max_shell_type)
        gb4w = new gbw.GB4IntegralWrapper(gobasis._this,
                            <ints.GB4Integral*> gb4int)
        gb4w.compute_diagonal(&output[0, 0])

    finally:
        if gb4int is not NULL:
            del gb4int
        if gb4w is not NULL:
            del gb4w

def _get_2index_slice(GOBasis gobasis, long index0, long index2,
                     double[:, ::1] index_slice not None):
    """Get a 2-index slice. For cholesky testing only."""
    cdef ints.GB4ElectronRepulsionIntegralLibInt* gb4int = NULL
    cdef gbw.GB4IntegralWrapper* gb4w = NULL
    assert index_slice.shape[0] == gobasis.nbasis
    assert index_slice.shape[1] == gobasis.nbasis

    cdef long pbegin0
    cdef long pend0
    cdef long pbegin2
    cdef long pend2
    cdef double* output
    try:
        gb4int = new ints.GB4ElectronRepulsionIntegralLibInt(
                            gobasis.max_shell_type)
        gb4w = new gbw.GB4IntegralWrapper(gobasis._this,
                            <ints.GB4Integral*> gb4int)
        gb4w.select_2index(index0, index2, &pbegin0, &pend0, &pbegin2, &pend2)
        gb4w.compute()
        output = gb4w.get_2index_slice(index0, index2)
        print(output[0])
        print(sizeof(double)*gobasis.nbasis*gobasis.nbasis)
        libc.string.memcpy(&index_slice[0, 0], output,
                           sizeof(double) * gobasis.nbasis * gobasis.nbasis)
        print(index_slice[0, 0])

    finally:
        if gb4int is not NULL:
            del gb4int
        if gb4w is not NULL:
            del gb4w

#
# ints wrappers (for testing only)
#

cdef class _GB2Integral:
    """Wrapper for ints.GB2Integral, for testing only"""
    cdef ints.GB2Integral* _baseptr

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

    def reset(self, long shell_type0, long shell_type1,
              double[::1] r0 not None,
              double[::1] r1 not None):
        assert r0.shape[0] == 3
        assert r1.shape[0] == 3
        self._baseptr.reset(shell_type0, shell_type1, &r0[0], &r1[0])

    def add(self, double coeff, double alpha0, double alpha1,
            double[::1] scales0 not None,
            double[::1] scales1 not None):
        assert scales0.shape[0] == _get_shell_nbasis(abs(self._baseptr.get_shell_type0()))
        assert scales1.shape[0] == _get_shell_nbasis(abs(self._baseptr.get_shell_type1()))
        self._baseptr.add(coeff, alpha0, alpha1, &scales0[0], &scales1[0])

    def cart_to_pure(self):
        self._baseptr.cart_to_pure()

    def get_work(self, shape0, shape1):
        """This returns a **copy** of the c++ work array.

           Returning a numpy array with a buffer created in c++ is dangerous.
           If the c++ array becomes deallocated, the numpy array may still
           point to the deallocated memory. For that reason, a copy is returned.
           Speed is not an issue as this class is only used for testing.
        """
        cdef np.npy_intp shape[2]
        assert shape0 > 0
        assert shape1 > 0
        assert shape0 <= self.max_nbasis
        assert shape1 <= self.max_nbasis
        shape[0] = shape0
        shape[1] = shape1
        tmp = np.PyArray_SimpleNewFromData(2, shape, np.NPY_DOUBLE, <void*> self._baseptr.get_work())
        return tmp.copy()


cdef class _GB2OverlapIntegral(_GB2Integral):
    """Wrapper for ints.GB2OverlapIntegral, for testing only"""
    cdef ints.GB2OverlapIntegral* _this

    def __cinit__(self, long max_nbasis):
        self._this = new ints.GB2OverlapIntegral(max_nbasis)
        self._baseptr = <ints.GB2Integral*> self._this


cdef class _GB2KineticIntegral(_GB2Integral):
    """Wrapper for ints.GB2KineticIntegral, for testing only"""
    cdef ints.GB2KineticIntegral* _this

    def __cinit__(self, long max_nbasis):
        self._this = new ints.GB2KineticIntegral(max_nbasis)
        self._baseptr = <ints.GB2Integral*> self._this


cdef class _GB2NuclearAttractionIntegral(_GB2Integral):
    """Wrapper for ints.GB2NuclearAttractionIntegral, for testing only"""
    # make an additional reference to these arguments to avoid deallocation
    cdef double[::1] _charges
    cdef double[:, ::1] _centers
    cdef ints.GB2NuclearAttractionIntegral* _this

    def __cinit__(self, long max_nbasis,
                  double[::1] charges not None,
                  double[:, ::1] centers not None):
        cdef long ncharge = charges.shape[0]
        assert centers.shape[0] == ncharge
        self._charges = charges
        self._centers = centers
        self._this = new ints.GB2NuclearAttractionIntegral(
           max_nbasis, &charges[0], &centers[0, 0], ncharge
        )
        self._baseptr = <ints.GB2Integral*> self._this


cdef class _GB2ErfAttractionIntegral(_GB2Integral):
    """Wrapper for ints.GB2ErfAttractionIntegral, for testing only"""
    # make an additional reference to these arguments to avoid deallocation
    cdef double[::1] _charges
    cdef double[:, ::1] _centers
    cdef ints.GB2ErfAttractionIntegral* _this

    def __cinit__(self, long max_nbasis,
                  double[::1] charges not None,
                  double[:, ::1] centers not None, double mu):
        cdef long ncharge = charges.shape[0]
        assert centers.shape[0] == ncharge
        self._charges = charges
        self._centers = centers
        self._this = new ints.GB2ErfAttractionIntegral(
            max_nbasis, &charges[0], &centers[0, 0], ncharge, mu
        )
        self._baseptr = <ints.GB2Integral*> self._this

    @property
    def mu(self):
        return self._this.get_mu()


cdef class _GB2GaussAttractionIntegral(_GB2Integral):
    """Wrapper for ints.GB2GaussAttractionIntegral, for testing only"""
    # make an additional reference to these arguments to avoid deallocation
    cdef double[::1] _charges
    cdef double[:, ::1] _centers
    cdef ints.GB2GaussAttractionIntegral* _this

    def __cinit__(self, long max_nbasis,
                  double[::1] charges not None,
                  double[:, ::1] centers not None, double c,
                  double alpha):
        cdef long ncharge = charges.shape[0]
        assert centers.shape[0] == ncharge
        self._charges = charges
        self._centers = centers
        self._this = new ints.GB2GaussAttractionIntegral(
            max_nbasis, &charges[0], &centers[0, 0], ncharge, c, alpha
        )
        self._baseptr = <ints.GB2Integral*> self._this

    @property
    def c(self):
        return self._this.get_c()

    @property
    def alpha(self):
        return self._this.get_alpha()


#
# iter_gb wrappers (for testing only)
#


cdef class _IterGB1:
    """Wrapper for the IterGB1 class, for testing only."""
    cdef iter_gb.IterGB1* _this
    cdef GBasis _gbasis

    def __cinit__(self, GBasis gbasis not None):
        self._this = new iter_gb.IterGB1(gbasis._baseptr)
        self._gbasis = gbasis

    def __dealloc__(self):
        del self._this

    def inc_shell(self):
        return self._this.inc_shell()

    def update_shell(self):
        self._this.update_shell()

    def inc_prim(self):
        return self._this.inc_prim()

    def update_prim(self):
        self._this.update_prim()

    def store(self, double[::1] work not None,
              double[::1] output not None, long dim=1):
        max_shell_nbasis = _get_shell_nbasis(self._gbasis.max_shell_type)
        assert work.shape[0] == _get_shell_nbasis(self._this.shell_type0)
        assert output.shape[0] == self._gbasis.nbasis
        self._this.store(&work[0], &output[0], dim)

    @property
    def public_fields(self):
        return (
            self._this.con_coeff,
            self._this.shell_type0,
            self._this.alpha0,
            self._this.r0[0], self._this.r0[1], self._this.r0[2],
            self._this.ibasis0,
        )

    @property
    def private_fields(self):
        return (
            self._this.ishell0,
            self._this.nprim0,
            self._this.oprim0,
            self._this.iprim0,
        )


cdef class _IterGB2:
    """Wrapper for the IterGB2 class, for testing only."""
    cdef iter_gb.IterGB2* _this
    cdef GBasis _gbasis

    def __cinit__(self, GBasis gbasis not None):
        self._this = new iter_gb.IterGB2(gbasis._baseptr)
        self._gbasis = gbasis

    def __dealloc__(self):
        del self._this

    def inc_shell(self):
        return self._this.inc_shell()

    def update_shell(self):
        self._this.update_shell()

    def inc_prim(self):
        return self._this.inc_prim()

    def update_prim(self):
        self._this.update_prim()

    def store(self, double[:, ::1] work not None,
              double[:, ::1] output not None):
        max_shell_nbasis = _get_shell_nbasis(self._gbasis.max_shell_type)
        assert work.shape[0] == _get_shell_nbasis(self._this.shell_type0)
        assert work.shape[1] == _get_shell_nbasis(self._this.shell_type1)
        assert output.shape[0] == self._gbasis.nbasis
        assert output.shape[1] == self._gbasis.nbasis
        self._this.store(&work[0, 0], &output[0, 0])

    @property
    def public_fields(self):
        return (
            self._this.con_coeff,
            self._this.shell_type0, self._this.shell_type1,
            self._this.alpha0, self._this.alpha1,
            self._this.r0[0], self._this.r0[1], self._this.r0[2],
            self._this.r1[0], self._this.r1[1], self._this.r1[2],
            self._this.ibasis0, self._this.ibasis1,
        )

    @property
    def private_fields(self):
        return (
            self._this.ishell0, self._this.ishell1,
            self._this.nprim0, self._this.nprim1,
            self._this.oprim0, self._this.oprim1,
            self._this.iprim0, self._this.iprim1,
        )


cdef class _IterGB4:
    """Wrapper for the IterGB4 class, for testing only."""
    cdef iter_gb.IterGB4* _this
    cdef GBasis _gbasis

    def __cinit__(self, GBasis gbasis not None):
        self._this = new iter_gb.IterGB4(gbasis._baseptr)
        self._gbasis = gbasis

    def __dealloc__(self):
        del self._this

    def inc_shell(self):
        return self._this.inc_shell()

    def update_shell(self):
        self._this.update_shell()

    def inc_prim(self):
        return self._this.inc_prim()

    def update_prim(self):
        self._this.update_prim()

    def store(self, double[:, :, :, ::1] work not None,
              double[:, :, :, ::1] output not None):
        max_shell_nbasis = _get_shell_nbasis(self._gbasis.max_shell_type)
        assert work.shape[0] == _get_shell_nbasis(self._this.shell_type0)
        assert work.shape[1] == _get_shell_nbasis(self._this.shell_type1)
        assert work.shape[2] == _get_shell_nbasis(self._this.shell_type2)
        assert work.shape[3] == _get_shell_nbasis(self._this.shell_type3)
        assert output.shape[0] == self._gbasis.nbasis
        assert output.shape[1] == self._gbasis.nbasis
        assert output.shape[2] == self._gbasis.nbasis
        assert output.shape[3] == self._gbasis.nbasis
        self._this.store(&work[0, 0, 0, 0], &output[0, 0, 0, 0])

    @property
    def public_fields(self):
        return (
            self._this.con_coeff,
            self._this.shell_type0, self._this.shell_type1, self._this.shell_type2, self._this.shell_type3,
            self._this.alpha0, self._this.alpha1, self._this.alpha2, self._this.alpha3,
            self._this.r0[0], self._this.r0[1], self._this.r0[2],
            self._this.r1[0], self._this.r1[1], self._this.r1[2],
            self._this.r2[0], self._this.r2[1], self._this.r2[2],
            self._this.r3[0], self._this.r3[1], self._this.r3[2],
            self._this.ibasis0, self._this.ibasis1, self._this.ibasis2, self._this.ibasis3,
        )

    @property
    def private_fields(self):
        return (
            self._this.ishell0, self._this.ishell1, self._this.ishell2, self._this.ishell3,
            self._this.nprim0, self._this.nprim1, self._this.nprim2, self._this.nprim3,
            self._this.oprim0, self._this.oprim1, self._this.oprim2, self._this.oprim3,
            self._this.iprim0, self._this.iprim1, self._this.iprim2, self._this.iprim3,
        )


#
# iter_pow wrappers (for testing only)
#


def _iter_pow1_inc(long[::1] n not None):
    assert n.shape[0] == 3
    return iter_pow.iter_pow1_inc(&n[0])


cdef class _IterPow1:
    """Wrapper for the IterPow1 class, for testing only."""
    cdef iter_pow.IterPow1* _c_i1p

    def __cinit__(self):
        self._c_i1p = new iter_pow.IterPow1()

    def __dealloc__(self):
        del self._c_i1p

    def __init__(self, long shell_type0):
        if shell_type0 < 0:
            raise ValueError('A shell_type parameter can not be negative.')
        self._c_i1p.reset(shell_type0)

    def inc(self):
        return self._c_i1p.inc()

    @property
    def fields(self):
        return (
            self._c_i1p.n0[0], self._c_i1p.n0[1], self._c_i1p.n0[2],
            self._c_i1p.ibasis0
        )


cdef class _IterPow2:
    """Wrapper for the IterPow2 class, for testing only."""
    cdef iter_pow.IterPow2* _c_i2p

    def __cinit__(self):
        self._c_i2p = new iter_pow.IterPow2()

    def __dealloc__(self):
        del self._c_i2p

    def __init__(self, long shell_type0, long shell_type1):
        if shell_type0 < 0 or shell_type1 < 0:
            raise ValueError('A shell_type parameter can not be negative.')
        self._c_i2p.reset(shell_type0, shell_type1)

    def inc(self):
        return self._c_i2p.inc()

    @property
    def fields(self):
        return (
            self._c_i2p.n0[0], self._c_i2p.n0[1], self._c_i2p.n0[2],
            self._c_i2p.n1[0], self._c_i2p.n1[1], self._c_i2p.n1[2],
            self._c_i2p.offset, self._c_i2p.ibasis0, self._c_i2p.ibasis1,
        )