#cython: language_level=3

from typing import Tuple, Union

cimport numpy as np
from gbasis.pxds cimport c_gbasis

# cdef _prepare_array(array: Union[np.ndarray, None], shape: Tuple[int], name: str)
# cdef _check_shape(array: Union[np.ndarray, memoryview], shape: Tuple[int], name: str)

cdef _prepare_array(array, shape, name)
cdef _check_shape(array, shape, name)
cpdef _get_shell_nbasis(long shell_type)

cdef class GBasis:
    cdef c_gbasis.GBasis* _baseptr
    # Keep reference to arrays to make sure they will not be deallocated.
    cdef np.ndarray _centers
    cdef np.ndarray _shell_map
    cdef np.ndarray _nprims
    cdef np.ndarray _shell_types
    cdef np.ndarray _alphas
    cdef np.ndarray _con_coeffs
