cimport gbasis
from cext_common cimport GBasis

cdef class GOBasis(GBasis):
    cdef public list _biblio
    cdef gbasis.GOBasis* _this