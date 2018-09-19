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
GBasis is a HORTON3 package for calculating Gaussian Basis integrals.

License
-------

GBasis is distributed under GPL License version 3 (GPLv3).


Dependencies
------------

Gbasis (like all HORTON3 packages) is built using conda:

It is provided in the conda **theochem** channel.


Installation (from Conda)
-------------------------

To install GBasis:

```bash
conda -c theochem install gbasis
```

Installation (from source)
--------------------------

If you wish to build from source, you will need the **conda-build** package
to build it.

You must set the PROJECT_VERSION and MYCONDAPY environmental variables to
emulate the travis build environment.

From project root, issue some variation of:

```bash
PROJECT_VERSION=0.0.0 MYCONDAPY=3.7 conda-build tools/conda.recipe/meta.yaml
```

Installation (manual)
---------------------

Advanced developers may build by hand using the dependencies listed below,
but the procedure is entirely unsupported. You are responsible for setting
the proper flags to find headers and linking.

The following dependencies will be necessary for GBasis to build properly,

* Python >= 3.6
* SciPy >= 0.11.0
* NumPy >= 1.9.1
* Nosetests >= 1.1.2
* Libint = 2.0.3
* gcc/clang with C++ 1x support


Testing
-------

The tests are automatically run when building with conda, but you may try
them again on your own machine:

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
