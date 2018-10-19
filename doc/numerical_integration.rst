Computing an integral involving the electron density
====================================================

This section assumes that the following objects are already available:

* ``obasis``: an orbital basis set object
* ``dm_full``: a spin-summed density matrix
* ``grid``: a Becke-Lebedev integration grid (available from the ``grids`` package).

If you are not familiar with the ``obasis`` and ``dm_full`` objects, please
refer to :ref:`user_molecularham_basis` and :ref:`user_hf_dft`,
respectively. Note that the density matrix can also be loaded from a file
instead of being computed by GBasis. This is demonstrated in an example at the
end of this section.

First, you must evaluate the integrand on all the points of the integration grid.
In case of the electron density, this can be done as follows:

.. code-block:: python

    rho = obasis.compute_grid_density_dm(dm_full, grid.points)

It is assumed that ``dm_full`` is symmetric. Several other quantities can also
be evaluated on the grid. For more details, please refer to:

* :py:meth:`~gbasis.cext.GOBasis.compute_grid_density_dm`
* :py:meth:`~gbasis.cext.GOBasis.compute_grid_gradient_dm`
* :py:meth:`~gbasis.cext.GOBasis.compute_grid_gga_dm`
* :py:meth:`~gbasis.cext.GOBasis.compute_grid_kinetic_dm`
* :py:meth:`~gbasis.cext.GOBasis.compute_grid_hessian_dm`
* :py:meth:`~gbasis.cext.GOBasis.compute_grid_mgga_dm`
* :py:meth:`~gbasis.cext.GOBasis.compute_grid_hartree_dm`
* :py:meth:`~gbasis.cext.GOBasis.compute_grid_esp_dm`
* :py:meth:`~gbasis.cext.GOBasis.compute_grid_orbitals_exp`

Integrating the electron density by itself results in the total number of electrons.
This is a simple way to verify the accuracy of the integration grid.

.. code-block:: python

    print(grid.integrate(rho))

Since ``rho`` is simply a Numpy array, it can be manipulated easily to compute
functionals of the electron density, e.g.

.. code-block:: python

    print(grid.integrate(rho**(4.0/3.0)))

You can also use the ``grid.points`` array to evaluate other expectation values
numerically, e.g. the following snippet evaluates the expectation value of
:math:`\vert\mathbf{r}\vert=(x^{2}+y^{2}+z^{2})^{0.5}`:

.. code-block:: python

    r = (grid.points[:,0]**2 + grid.points[:,1]**2 + grid.points[:,2]**2)**0.5
    grid.integrate(rho, r)

As shown in the above snippet, the ``integrate`` method can take multiple
one-dimensional arrays that are all multiplied before integration.

The following script is a complete example for computing the expectation value
of :math:`\vert\mathbf{r}\vert=(x^{2}+y^{2}+z^{2})^{0.5}`
for a molecular wave-function loaded from a file.

.. literalinclude:: ../data/examples/grid/expectation_r.py
    :caption: data/examples/grid/expectation_r.py


Constructing a one-body operator from a real-space potential
============================================================

This section assumes that the following objects are already available:

* ``obasis``: an orbital basis set object
* ``dm_full``: a spin-summed density matrix
* ``grid``: a Becke-Lebedev integration grid as introduced above.

If you are not familiar with the ``obasis`` object, please refer to
:ref:`user_molecularham_basis`. The density matrix can either be read from a
file or computed with GBasis. For more information, please refer to :ref:`user_hf_dft`.

Given a multiplicative potential, its expectation value is written as:

.. math::

    \langle V \rangle = \int \rho(\mathbf{r}) V(\mathbf{r}) d\mathbf{r}.

Expanding the orbitals in a local basis set results in:

.. math::

    \langle V \rangle = \sum_{\mu\nu} D_{\mu\nu} \mathcal{V}_{\nu\mu}

where :math:`D_{\mu\nu}` is the spin-summed density matrix. The matrix
:math:`\mathcal{V}_{\nu\mu}` is defined as

.. math::

    \mathcal{V}_{\nu\mu} = \int V(\mathbf{r}) b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) d\mathbf{r}

where :math:`b_\mu(\mathbf{r})` are the orbital basis functions. Such matrices
can be constructed with the
:py:meth:`~gbasis.cext.GOBasis.compute_grid_density_fock` method. This method is
also useful when applying the chain rule to construct the contribution of a density functional
to a Fock matrix:

.. math::

    \frac{\partial E[\rho]}{\partial D_{\nu\mu}} = \int \frac{\delta E[\rho]}{\delta \rho(\mathbf{r})}  b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) d\mathbf{r}

The usage pattern is as follows:

.. code-block:: python

    # Construct some potential, e.g. a hyperbolic well
    rsq = grid.points[:,0]**2 + grid.points[:,1]**2 + grid.points[:,2]**2
    pot = np.sqrt(1 + rsq)

    # Actual computation
    fock = obasis.compute_grid_density_fock(grid.points, grid.weights, pot)


Other chain rules are also implemented:

* :py:meth:`~gbasis.cext.GOBasis.compute_grid_gradient_fock`

  .. math::

    \frac{\partial E[\nabla\rho]}{\partial D_{\nu\mu}} =
        \int \frac{\delta E[\nabla\rho]}{\delta \nabla\rho(\mathbf{r})}
        \cdot \left(
            \nabla b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) +
            b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r})
        \right) d\mathbf{r}

* :py:meth:`~gbasis.cext.GOBasis.compute_grid_gradient_fock` just
  combines the density and gradient chain rules. This is more efficient than
  computing the separately.

* :py:meth:`~gbasis.cext.GOBasis.compute_grid_kinetic_fock`

  .. math::

    \frac{\partial E[\tau]}{\partial D_{\nu\mu}} =
        \frac{1}{2}\int \frac{\delta E[\tau]}{\delta \tau(\mathbf{r})}
        \nabla b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r}) d\mathbf{r}

  where :math:`\tau(\mathbf{r})` is the positive kinetic energy density:

  ..
    Mind the adjective "positive" in the following sentence. There are many
    choices for the kinetic energy density (which is essentially arbitrarily
    defined). Technically one should say "the nonnegative kinetic energy
    density" but I think most people call this the "positive" choice.
        ~ Paul W. Ayers

  .. math::

    \tau(\mathbf{r}) = \frac{1}{2} \sum_{\mu\nu} D_{\mu\nu} \nabla b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r})

* :py:meth:`~gbasis.cext.GOBasis.compute_grid_hessian_fock`

  .. math::

    \frac{\partial E[\nabla\nabla\rho]}{\partial D_{\nu\mu}} =
        \int \frac{\delta E[\nabla\nabla\rho]}{\delta \nabla\nabla\rho(\mathbf{r})}
        \colon \left(
            \nabla \nabla b_\nu^*(\mathbf{r}) b_\mu(\mathbf{r}) +
            2 \nabla b_\nu^*(\mathbf{r}) \nabla b_\mu(\mathbf{r}) +
            b_\nu^*(\mathbf{r}) \nabla \nabla b_\mu(\mathbf{r})
        \right) d\mathbf{r}

* :py:meth:`~gbasis.cext.GOBasis.compute_grid_mgga_fock` just combines
  several of the previous chain rules: density, gradient, laplacian (trace of
  the Hessian) and kinetic energy density.