#!/usr/bin/env bash

if [[ -z "${BUILD_DEBUG}" ]]; then
    mkdir build_release
    cd build_release
    cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ..
    make && make install
else
    mkdir build_debug
    cd build_debug
    cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_BUILD_TYPE=debug ..
    VERBOSE=1 make && make install
fi
