import numpy as np

from .common import (load_obasis1, load_obasis2, load_obasis_grid, load_mdata, load_dm, load_er,
                     load_quad, load_dipole, load_na, load_kin, load_olp)
from .lightgrid import generate_molecular_grid, integrate


def check_g09_overlap(fn):
    obasis = load_obasis1(fn)
    olp1 = obasis.compute_overlap()
    olp2 = load_olp(fn)
    mask = abs(olp1) > 1e-5
    delta = olp1 - olp2
    expect = olp1
    error = (delta[mask] / expect[mask]).max()
    assert error < 1e-5


def test_overlap_water_sto3g_hf():
    check_g09_overlap('water_sto3g_hf_g03_fchk')


def test_overlap_water_ccpvdz_pure_hf():
    check_g09_overlap('water_ccpvdz_pure_hf_g03_fchk')


def test_overlap_water_ccpvdz_cart_hf():
    check_g09_overlap('water_ccpvdz_cart_hf_g03_fchk')


def test_overlap_co_ccpv5z_pure_hf():
    check_g09_overlap('co_ccpv5z_pure_hf_g03_fchk')


def test_overlap_co_ccpv5z_cart_hf():
    check_g09_overlap('co_ccpv5z_cart_hf_g03_fchk')


def check_g09_kinetic(fn):
    obasis = load_obasis1(fn)
    kin1 = obasis.compute_kinetic()
    kin2 = load_kin(fn)
    mask = abs(kin1) > 1e-5
    delta = kin1 - kin2
    expect = kin1
    error = (delta[mask] / expect[mask]).max()
    assert error < 1e-5


def test_kinetic_water_sto3g_hf():
    check_g09_kinetic('water_sto3g_hf_g03_fchk')


def test_kinetic_water_ccpvdz_pure_hf():
    check_g09_kinetic('water_ccpvdz_pure_hf_g03_fchk')


def test_kinetic_water_ccpvdz_cart_hf():
    check_g09_kinetic('water_ccpvdz_cart_hf_g03_fchk')


def test_kinetic_co_ccpv5z_pure_hf():
    check_g09_kinetic('co_ccpv5z_pure_hf_g03_fchk')


def test_kinetic_co_ccpv5z_cart_hf():
    check_g09_kinetic('co_ccpv5z_cart_hf_g03_fchk')


def check_g09_nuclear_attraction(fn):
    mol = load_mdata(fn)
    obasis = load_obasis1(fn)
    na1 = obasis.compute_nuclear_attraction(mol['coordinates'], mol['pseudo_numbers'])
    na2 = load_na(fn)
    mask = abs(na1) > 1e-5
    expect = na1
    result = na2
    delta = -expect - result
    error = (delta[mask] / expect[mask]).max()
    assert error < 4e-5


def test_nuclear_attraction_water_sto3g_hf():
    check_g09_nuclear_attraction('water_sto3g_hf_g03_fchk')


def test_nuclear_attraction_water_ccpvdz_pure_hf():
    check_g09_nuclear_attraction('water_ccpvdz_pure_hf_g03_fchk')


def test_nuclear_attraction_water_ccpvdz_cart_hf():
    check_g09_nuclear_attraction('water_ccpvdz_cart_hf_g03_fchk')


def test_nuclear_attraction_co_ccpv5z_pure_hf():
    check_g09_nuclear_attraction('co_ccpv5z_pure_hf_g03_fchk')


def test_nuclear_attraction_co_ccpv5z_cart_hf():
    check_g09_nuclear_attraction('co_ccpv5z_cart_hf_g03_fchk')


def check_g09_dipole(fn):
    """Compare dipole moment computed from WFN and nuclei to reference value.

    Parameters
    ----------
    fn : str
        The sanitized FCHK filename with '.' and '-' substituted with '_'.

    """
    obasis = load_obasis1(fn)
    mol = load_mdata(fn)
    xyz_array = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    center = np.zeros(3)
    dm_full = load_dm(fn)
    dipole = []
    for xyz in xyz_array:
        dipole_ints = obasis.compute_multipole_moment(xyz, center)
        dipole_v = -np.einsum('ab,ba', dipole_ints, dm_full)
        for i, n in enumerate(mol['pseudo_numbers']):
            dipole_v += n * pow(mol['coordinates'][i, 0], xyz[0]) * \
                        pow(mol['coordinates'][i, 1], xyz[1]) * \
                        pow(mol['coordinates'][i, 2], xyz[2])
        dipole.append(dipole_v)
    np.testing.assert_almost_equal(dipole, load_dipole(fn), decimal=6)


def test_dipole_water_sto3g_hf():
    check_g09_dipole('water_sto3g_hf_g03_fchk')


def test_dipole_water_ccpvdz_pure_hf():
    check_g09_dipole('water_ccpvdz_pure_hf_g03_fchk')


def test_dipole_ccpvdz_cart_hf():
    check_g09_dipole('water_ccpvdz_cart_hf_g03_fchk')


def test_dipole_co_ccpv5z_pure_hf():
    check_g09_dipole('co_ccpv5z_pure_hf_g03_fchk')


def test_dipole_co_ccpv5z_cart_hf():
    check_g09_dipole('co_ccpv5z_cart_hf_g03_fchk')


def check_g09_quadrupole(fn):
    """Compare quadrupole moment computed from WFN and nuclei to reference value.

    Parameters
    ----------
    fn : str
        The sanitized FCHK filename with '.' and '-' substituted with '_'.

    """
    obasis = load_obasis1(fn)
    mol = load_mdata(fn)
    # HORTON ordering: xx, xy, xz, yy, yz, zz (alphabetic)
    xyz_array = np.array([[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]])
    center = np.zeros(3)
    dm_full = load_dm(fn)
    quadrupole = []
    for xyz in xyz_array:
        quadrupole_ints = obasis.compute_multipole_moment(xyz, center)
        quad_v = -np.einsum('ab,ba', quadrupole_ints, dm_full)
        for i, n in enumerate(mol['pseudo_numbers']):
            quad_v += n * pow(mol['coordinates'][i, 0], xyz[0]) * \
                      pow(mol['coordinates'][i, 1], xyz[1]) * \
                      pow(mol['coordinates'][i, 2], xyz[2])
        quadrupole.append(quad_v)
    # removing trace:
    mean = (quadrupole[0] + quadrupole[3] + quadrupole[5]) / 3.0
    quadrupole[0] -= mean
    quadrupole[3] -= mean
    quadrupole[5] -= mean
    np.testing.assert_almost_equal(quadrupole, load_quad(fn), decimal=6)


def test_quadrupole_ch3_hf_sto3g():
    check_g09_quadrupole('ch3_hf_sto3g_fchk')


def test_quadrupole_li_h_321g_hf_g09():
    check_g09_quadrupole('li_h_3_21G_hf_g09_fchk')


def test_quadrupole_monosilicic_acid_hf_lan_g09():
    check_g09_quadrupole('monosilicic_acid_hf_lan_fchk')


def check_g09_electron_repulsion(fn, check_g09_zeros=False):
    obasis = load_obasis2(fn)
    er1 = obasis.compute_electron_repulsion()
    er2 = load_er(fn)
    mask = abs(er1) > 1e-6
    expect = er1
    got = er2
    if check_g09_zeros:
        assert ((expect == 0.0) == (got == 0.0)).all()
    delta = expect - got
    error = (delta[mask] / expect[mask]).max()
    assert error < 1e-5


def test_electron_repulsion_water_sto3g_hf():
    check_g09_electron_repulsion('water_sto3g_hf_g03_fchk', True)


def test_electron_repulsion_water_ccpvdz_pure_hf():
    check_g09_electron_repulsion('water_ccpvdz_pure_hf_g03_fchk')


def test_electron_repulsion_water_ccpvdz_cart_hf():
    check_g09_electron_repulsion('water_ccpvdz_cart_hf_g03_fchk')


def check_g09_grid_fn(fn):
    obasis = load_obasis_grid(fn)
    obasis1 = load_obasis1(fn)
    mol = load_mdata(fn)
    points, weights = generate_molecular_grid(mol['numbers'], mol['coordinates'], 10000)
    dm_full = load_dm(fn)
    rhos = obasis.compute_grid_density_dm(dm_full, points)
    pop = integrate(weights, rhos)
    nel = np.einsum('ab,ba', obasis1.compute_overlap(), dm_full)
    np.testing.assert_allclose(pop, nel, atol=1e-2)


def test_grid_fn_h_sto3g():
    check_g09_grid_fn('h_sto3g_fchk')


def test_grid_fn_lih_321g_hf():
    check_g09_grid_fn('li_h_3_21G_hf_g09_fchk')


def test_grid_fn_water_sto3g_hf_t():
    check_g09_grid_fn('water_sto3g_hf_g03_fchk')


def test_grid_fn_co_ccpv5z_pure_hf_t():
    check_g09_grid_fn('co_ccpv5z_pure_hf_g03_fchk')


def test_grid_fn_co_ccpv5z_cart_hf_t():
    check_g09_grid_fn('co_ccpv5z_cart_hf_g03_fchk')
