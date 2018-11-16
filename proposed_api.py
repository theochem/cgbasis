# Proposed GBasis API after splitting C++/Python, Ones, Twos, Grid

# File structure (assume standard docs, conda, setup.py, cmake are present...):
# gbasis/
# - gbasis.[h/cpp]
# - common.[h/cpp]
# - calc.[h/cpp]
# - cartpure.[h/cpp]
# - ones/
# - - iter_pow1.[h/cpp]
# - - iter_gb1.[h/cpp]
# - - nucpot.[h/cpp]
# - - ints1.[h/cpp]
# - twos/
# - - iter_gb2.[h/cpp]
# - - ints2.[h/cpp]
# - - boys.[h/cpp]
# - - gbw.[h/cpp]
# - grid/
# - - fns.[h/cpp]
# - - iter_pow_grid.[h/cpp]
# - sparse/
# - - cholesky.[h/cpp]
# python/
# - test/
# - - (all test ext)
# - - (all test pys)
# - ext/
# - - (all non-test ext, ie c_gbasis.pxd)
# - periodic.py
# - utils.py
# - gobasis.py
# - ones.pyx
# - twos.pyx
# - grid.pyx
# - sparse.pyx

# The API would then become:
gbasis_params = {...}
gb1 = gbasis.get_gbasis1(**gbasis_params)
gb2 = gbasis.get_gbasis2(**gbasis_params)
gbg = gbasis.get_gbasisgrid(**gbasis_params)
gbs = gbasis.get_gbasis_sparse(**gbasis_params)

olp = gb1.compute_overlap()
kin = gb1.compute_overlap()
# ...

eri = gb2.compute_electron_repulsion()
eri_alpha = gb2.compute_ralpha_repulsion()

grid = gbg.compute_grid_density_dm(...)

cholesky2e = gbs.compute_cholesky_electron_repulsion()