GBasis
======
|Travis|
|Conda|
|Pypi|
|Codecov|
|Version|
|CondaVersion|

About
-----


License
-------

GBasis is distributed under GPL License version 3 (GPLv3).


Dependencies
------------

The following dependencies will be necessary for GBasis to build properly,

* Python >= 2.7, >= 3.5: http://www.python.org/
* SciPy >= 0.11.0: http://www.scipy.org/
* NumPy >= 1.9.1: http://www.numpy.org/
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/


Installation
------------

To install GBasis:

```bash
./tools/libs/install_libint-2.0.3.sh
python ./setup install --user
```


Testing
-------

To run tests:

```bash
nosetests -v gbasis
```

.. |Travis| image:: https://travis-ci.org/theochem/gbasis.svg?branch=master
    :target: https://travis-ci.org/theochem/gbasis
.. |Version| image:: https://img.shields.io/pypi/pyversions/gbasis.svg
.. |Pypi| image:: https://img.shields.io/pypi/v/gbasis.svg
    :target: https://pypi.python.org/pypi/gbasis/0.1.3
.. |Codecov| image:: https://img.shields.io/codecov/c/github/theochem/gbasis/master.svg
    :target: https://codecov.io/gh/theochem/gbasis
.. |Conda| image:: https://img.shields.io/conda/v/theochem/gbasis.svg
    :target: https://anaconda.org/theochem/gbasis
.. |CondaVersion| image:: https://img.shields.io/conda/pn/theochem/gbasis.svg
    :target: https://anaconda.org/theochem/gbasis
