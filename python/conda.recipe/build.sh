${PYTHON} setup.py build_ext -I ${PREFIX}/include/gbasis:${PREFIX}/include/libint2 -L ${PREFIX}/lib
#PYTHONPATH=$PYTHONPATH:$PREFIX/lib/python3.7/site-packages python setup.py install --single-version-externally-managed --record=record.txt --prefix="${PREFIX}"
${PYTHON} setup.py install --prefix="${PREFIX}"
