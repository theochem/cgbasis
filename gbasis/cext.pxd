cimport numpy as np
cimport gbasis

cpdef _get_shell_nbasis(long shell_type)
cdef class GBasis:
    cdef gbasis.GBasis* _baseptr
    # Keep reference to arrays to make sure they will not be deallocated.
    cdef np.ndarray _centers
    cdef np.ndarray _shell_center
    cdef np.ndarray _shell_nprims
    cdef np.ndarray _shell_types
    cdef np.ndarray _alphas
    cdef np.ndarray _con_coeffs
cdef class GOBasis(GBasis):
    cdef public list _biblio
    cdef gbasis.GOBasis* _this