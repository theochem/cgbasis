Getting Started
===============

Using gbasis is fairly easy. The most common method is to use the
:code:`get_gobasis` method. Provide the coordinates, atomic numbers, and
basis set and it will return an instance of :code:`GOBasis`, which can be
used to get integrals.

.. code-block:: python

    import numpy as np
    from gbasis import get_gobasis

    coords = np.array([[0,0,0],
                       [0,0,1.47]])
    atnums = np.array([1,1])
    basis = "sto-3g"

    gb = get_gobasis(coords, atnums, basis)
    gb.compute_overlap()

.. warning::
    All HORTON packages use Bohr for coordinates.

If the basis set specified is not packaged, you may provide your own
in a file:

.. code-block:: python

    from gbasis import GOBasisFamily

    fam = GOBasisFamily("MyBasisSet*", "path/to/basis.nwchem")
    gb = get_gobasis(coords, atnums, fam)

The code can accept ``nwchem`` or ``gbs`` basis sets.

The library can also specify basis sets on particular atoms.
See :ref:`python-api` for more usage scenarios.

Calculating Integrals
---------------------

``Gbasis`` provides a variety of one and two electron integrals.

One electron integrals (coded in C++):

* Overlap
* Kinetic
* Nuclear Attraction
* Erf Attraction
* Gauss Attraction
* Multipole Moment

Two electron integrals (using LibInt):

* Electron Repulsion
* Erf Repulsion
* Gauss Repulsion
* RAlpha Repulsion
* Cholesky Electron Repulsion
* Cholesky Erf Repulsion
* Cholesky Gauss Repulsion
* Cholesky RAlpha Repulsion

Functions on grids
------------------

Additionally, ``GBasis`` can compute various functions on a
user-defined grid.

* Orbitals
* Orbital gradients
* Electron density (via density matrix)
* Electron density gradient (via density matrix)
* Positive definite kinetic energy density (via density matrix)
* Electron density Hessian (via density matrix)
* Electron density Laplacian (via density matrix)
* Hartree potential (via density matrix)
* Electrostatic potential (via density matrix)
*


