#cython: language_level=3

from gbasis.pxds cimport c_gbasis
from gbasis.cext_common cimport GBasis

cdef class GOBasis1(GBasis):
    cdef public list _biblio
    cdef c_gbasis.GOBasis* _this
