#!/usr/bin/env bash

mkdir build_release
cd build_release
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} ..
make && make install

