HORTON 3 Changelog
==================

The following API incompatible changes have been made:

* The functions only used for tests have been renamed private "_" in the gobasis.cext module.
* The Cython baseclasses no longer have the :code:`_this` pointer to the C++ object, instead it
  uses the :code:`_baseptr` pointer. This is to avoid the use of typecasts everywhere.
* 