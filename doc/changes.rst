HORTON 3 Changelog
==================

The following API incompatible changes have been made:

* The functions only used for tests have been renamed private "_" in the gobasis.cext module.
* The Cython baseclasses no longer have the :code:`_this` pointer to the C++ object, instead it
  uses the :code:`_baseptr` pointer. This is to avoid the use of typecasts everywhere.
* GOBasis has been split into 4 objects: GOBasis1, GOBasis2, GOBasisGrid, and GOBasisSparse.
  Each class now only contains the objects which pertain to it (i.e. 1e methods in GOBasis1).


The following *internal* changes have been made:

* The C++ library is now built (using cmake) separately from the Cython/Python library. Testing is
  still performed using Cython wrapped functions that end up using nosetests.
* The :code:`_compute_grid_point1` function has been moved from the GBasis class into the GOBasisGrid
  class.