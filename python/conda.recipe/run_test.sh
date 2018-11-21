#!/usr/bin/env bash
set -x

conda inspect linkages python-gbasis

pip install --upgrade git+https://github.com/theochem/cardboardlint.git@master#egg=cardboardlint
pip install pylint codecov coverage pycodestyle pydocstyle

if [[ -z "${BUILD_DEBUG}" ]]; then
    echo "Running Release tests"
    nosetests -v gbasis --detailed-errors
else
    echo "Running Debug tests"
    nosetests gbasis -v --detailed-errors --with-coverage --cover-package=gbasis --cover-tests --cover-branches
    coverage xml -i
fi