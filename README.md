GBasis
======
[![Build Status](https://travis-ci.org/theochem/gbasis.svg?branch=master)](https://travis-ci.org/theochem/gbasis)
<a href='https://docs.python.org/2.7/'><img src='https://img.shields.io/badge/python-2.7-blue.svg'></a>
<a href='https://docs.python.org/3.5/'><img src='https://img.shields.io/badge/python-3.5-blue.svg'></a>


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
