#!/usr/bin/env python

from iodata import IOData
from old_grids import BeckeMolGrid

from gbasis import get_gobasis

# Load the Gaussian output from file from HORTON's test data directory.
fn_fchk = 'water_sto3g_hf_g03.fchk'

# Replace the previous line with any other fchk file, e.g. fn_fchk = 'yourfile.fchk'.
mol = IOData.from_file(fn_fchk)

# Generate obasis instance from IOData
obasis = get_gobasis(**mol.obasis)

# Specify the integration grid
grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers)

# Get the spin-summed density matrix
dm_full = mol.get_dm_full()

# Compute the density on the grid
rho = obasis.compute_grid_density_dm(dm_full, grid.points)

# Compute the expectation value of |r|.
r = (grid.points[:, 0]**2 + grid.points[:, 1]**2 + grid.points[:, 2]**2)**0.5
expt_r = grid.integrate(rho, r)
print(f'EXPECTATION VALUE OF |R|: {expt_r}')
