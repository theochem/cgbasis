Install
=======

Dependencies
------------

Gbasis (like all HORTON3 packages) is built using conda:

It is provided in the conda **theochem** channel.


Installation (from Conda)
-------------------------

To install GBasis:

.. code-block:: bash

    $ conda -c theochem install gbasis

.. _install_from_source:

Installation (from source)
--------------------------

If you wish to build from source, you will need the **conda-build** package
to build it.

You must set the PROJECT_VERSION and MYCONDAPY environmental variables to
emulate the travis build environment.

From project root, issue some variation of:

.. code-block:: bash

    $ PROJECT_VERSION=0.0.0 MYCONDAPY=3.7 conda-build -c theochem tools/conda.recipe

Installation (by-hand)
----------------------

Advanced developers may build by hand using the dependencies listed below,
but the procedure is entirely unsupported. You are responsible for setting
the proper flags to find headers and linking through the setuptools -I, -L cli options.

The following dependencies will be necessary for GBasis to build properly,

* Python >= 3.7
* SciPy >= 0.11.0
* NumPy >= 1.9.1
* Nosetests
* Libint = 2.0.3
* gcc/clang with C++ 1x support
* Cython


Testing
-------

The tests are automatically run when building with conda, but you may try
them again on your own machine:

.. code-block:: bash

    $ nosetests -v gbasis
