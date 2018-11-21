#!/usr/bin/env bash

if [[ -z "${BUILD_DEBUG}" ]]; then
    ${PYTHON} setup.py build_ext -I ${PREFIX}/include/gbasis:${PREFIX}/include/libint2 \
        -L ${PREFIX}/lib
else
    ${PYTHON} setup.py build_ext -v --define CYTHON_TRACE_NOGIL \
        -I ${PREFIX}/include/gbasis:${PREFIX}/include/libint2 -L ${PREFIX}/lib
fi

${PYTHON} setup.py install --prefix="${PREFIX}"
