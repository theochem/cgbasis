from gbasis cimport iter_pow
from gbasis cimport iter_gb
from gbasis.cext cimport GBasis, _get_shell_nbasis

__all__ = [
    # iter_gb
    '_IterGB1', '_IterGB2', '_IterGB4',
    # iter_pow
    '_iter_pow1_inc', '_IterPow1', '_IterPow2',
]

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