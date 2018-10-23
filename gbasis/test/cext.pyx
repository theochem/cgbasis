from gbasis cimport iter_pow

__all__ = [
    # iter_pow
    '_iter_pow1_inc', '_IterPow1', '_IterPow2',
]


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