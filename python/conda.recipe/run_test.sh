#!/usr/bin/env bash
set -x

conda inspect linkages python-gbasis

if [[ -z "${BUILD_DEBUG}" ]]; then
    echo "Running Release tests"
    nosetests -v gbasis --detailed-errors
else
    echo "Running Debug tests"
    pip install --upgrade codecov coverage
    nosetests gbasis -v --detailed-errors --with-coverage --cover-package=gbasis --cover-tests --cover-branches
    coverage xml -i
fi