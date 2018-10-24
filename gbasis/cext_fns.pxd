cimport fns

cdef class _GB1DMGridFn:
    cdef fns.GB1DMGridFn* _baseptr

cdef class _GB1DMGridDensityFn(_GB1DMGridFn):
    cdef fns.GB1DMGridDensityFn* _this

cdef class _GB1DMGridGradientFn(_GB1DMGridFn):
    cdef fns.GB1DMGridGradientFn* _this

cdef class _GB1DMGridGGAFn(_GB1DMGridFn):
    cdef fns.GB1DMGridGGAFn* _this

cdef class _GB1DMGridKineticFn(_GB1DMGridFn):
    cdef fns.GB1DMGridKineticFn* _this

cdef class _GB1DMGridHessianFn(_GB1DMGridFn):
    cdef fns.GB1DMGridHessianFn* _this

cdef class _GB1DMGridMGGAFn(_GB1DMGridFn):
    cdef fns.GB1DMGridMGGAFn* _this