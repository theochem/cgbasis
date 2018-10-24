cimport ints

cdef class _GB4Integral:
    cdef ints.GB4Integral* _baseptr

cdef class _GB4ElectronRepulsionIntegralLibInt(_GB4Integral):
    cdef ints.GB4ElectronRepulsionIntegralLibInt* _this

cdef class _GB4ErfIntegralLibInt(_GB4Integral):
    cdef ints.GB4ErfIntegralLibInt* _this

cdef class _GB4GaussIntegralLibInt(_GB4Integral):
    cdef ints.GB4GaussIntegralLibInt* _this

cdef class _GB4RAlphaIntegralLibInt(_GB4Integral):
    cdef ints.GB4RAlphaIntegralLibInt* _this

cdef class _GB4DeltaIntegralLibInt(_GB4Integral):
    cdef ints.GB4DeltaIntegralLibInt* _this

cdef class _GB4IntraDensIntegralLibInt(_GB4Integral):
    cdef ints.GB4IntraDensIntegralLibInt* _this