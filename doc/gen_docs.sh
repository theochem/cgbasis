#!/bin/bash

PYTHONPATH="${PYTHONPATH}":.. python update_basis.py

doxygen

SPHINX_APIDOC_OPTIONS=members,undoc-members,show-inheritance,inherited-members sphinx-apidoc -o pyapi/ ../gbasis/ ../gbasis/test/test_*.py ../gbasis/bsets ../gbasis/test/cached --separate
