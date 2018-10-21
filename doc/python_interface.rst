.. _python-interface:

Python Interface
================

The user Python interface is fairly minimal, with the majority of the code
implemented in C++. There are substantial bits of Python code for testing
the C++ machinery, but they have been marked as private. The full user Python API
is documented in :ref:`pyapi`. If you would like to extend the Python API, see :ref:`developer`.

Some common operations using the Python API:

Unique basis set for all atoms
------------------------------

Usually, you want to define a unique (the same) basis set for the whole system.
This can be done by a function call.

.. code-block:: python

    obasis = get_gobasis(coordinates, numbers, 'cc-pvdz')

where ``mol.coordinates`` and ``mol.numbers`` are numpy arrays (see
:py:class:`~gbasis.cext.GOBasis`), and ``cc-pvdz`` is the cc-pVDZ basis set.

HORTON is distributed with most of the popular basis sets. A list of currently
supported built-in basis sets can be found here:
:ref:`ref_gaussian_basis_standard_sets`. The basis set for a given molecule is
constructed with the function :py:func:`~gbasis.gobasis.get_gobasis`


Specifying different basis sets for different atoms
---------------------------------------------------

In some cases, you may want to specify different basis sets for different atoms. For
example, you might like to use the 3-21G basis set for the hydrogen atom, the
6-31G basis set for the carbon atom, and STO-3G for all remaining atoms:

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',
                         element_map={'H':'3-21g', 'C':'6-31g'})

where `mol.coordinates` and `mol.numbers` are read from file (see
:py:class:`~gbasis.cext.GOBasis`),  and ``sto-3g``, ``3-21g`` and ``6-31g`` are the basis
set names (see ::ref:`ref_gaussian_basis_standard_sets`).

Alternatively, the same result can be obtained by substituting the H and C symbols
with their atomic numbers:

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',
                         element_map={1:'3-21g', 6:'6-31g'})

You can also override the default basis for selected atoms based on their index,
i.e. position in the list of atoms that specify the molecule:

.. code-block:: python

    obasis = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g',
                         index_map={0:'3-21g', 2:'6-31g'})

The above example uses the ``3-21g`` basis for the first atom, the ``6-31g``
basis for the third atom and the ``sto-3g`` basis for all other atoms.


Loading custom basis sets from file
-----------------------------------

You can also use other basis sets besides the ones that are shipped with HORTON.
It is assumed that the basis is available in NWChem format:

.. code-block:: python

    mybasis = GOBasisFamily('myname', filename='mybasis.nwchem'),
    obasis = get_gobasis(mol.coordinates, mol.numbers, mybasis)

Anywhere you can specify a built-in basis set with a string, you can also use
instance of the ``GOBasisFamily`` class (``mybasis`` in the example above), e.g.
in the arguments ``default``, ``element_map`` and ``index_map`` of
``get_gobasis``.


Defining basis sets with Python code
------------------------------------

In some circumstances, it may be useful to generate the basis set with some
Python code. For example, the following code generates an even tempered basis
for Lithium (without polarization functions):

.. literalinclude:: examples/even_tempered_li.py
    :caption: data/examples/hamiltonian/even_tempered_li.py

All basis functions in this example are just single s-type primitives, i.e. no
contactions are used. At the end of the example, the basis set is constructed
for a single Li atom in the origin.

Note that ``get_gobasis`` also accepts instances of GOBasisAtom for the
arguments ``default``, ``element_map`` and ``index_map``.

