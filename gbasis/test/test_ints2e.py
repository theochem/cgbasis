# -*- coding: utf-8 -*-
# pylint: skip-file
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


import numpy as np

from gbasis.cext import _get_shell_nbasis
from .common import (load_obasis, load_json)
from .. import (_GB4RAlphaIntegralLibInt, _GB4ErfIntegralLibInt, _GB4GaussIntegralLibInt,
                _GB4ElectronRepulsionIntegralLibInt,
                _GB4DeltaIntegralLibInt,
                _GB4IntraDensIntegralLibInt, get_gobasis)


def test_gb4_erilibint_class():
    max_shell_type = 4
    max_nbasis = _get_shell_nbasis(max_shell_type)

    gb4i = _GB4ElectronRepulsionIntegralLibInt(max_shell_type)
    assert gb4i.max_shell_type == max_shell_type
    assert gb4i.max_nbasis == max_nbasis
    assert gb4i.nwork == max_nbasis ** 4


def check_electron_repulsion(alphas0, alphas1, alphas2, alphas3, r0, r1, r2, r3, scales0, scales1,
                             scales2, scales3, shell_type0, shell_type1, shell_type2, shell_type3,
                             result0):
    # This test compares output from HORTON with reference data computed with
    # PyQuante.
    max_shell_type = 4
    max_nbasis = _get_shell_nbasis(max_shell_type)
    gb4i = _GB4ElectronRepulsionIntegralLibInt(max_shell_type)
    assert gb4i.max_nbasis == max_nbasis
    assert gb4i.nwork == max_nbasis ** 4

    nbasis0 = _get_shell_nbasis(shell_type0)
    nbasis1 = _get_shell_nbasis(shell_type1)
    nbasis2 = _get_shell_nbasis(shell_type2)
    nbasis3 = _get_shell_nbasis(shell_type3)
    assert result0.shape == (nbasis0, nbasis1, nbasis2, nbasis3)
    # Clear the working memory
    gb4i.reset(shell_type0, shell_type1, shell_type2, shell_type3, r0, r1, r2, r3)
    # Add a few contributions:
    for alpha0, alpha1, alpha2, alpha3 in zip(alphas0, alphas1, alphas2, alphas3):
        gb4i.add(1.0, alpha0, alpha1, alpha2, alpha3, scales0, scales1, scales2, scales3)
    result1 = gb4i.get_work(nbasis0, nbasis1, nbasis2, nbasis3)
    assert abs(result1 - result0).max() < 3e-7


def test_electron_repulsion_0_0_0_0_simple0():
    check_electron_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[4.37335457]]]]))


def test_electron_repulsion_0_0_0_0_simple1():
    check_electron_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[2.20567322]]]]))


def test_electron_repulsion_0_0_0_0_simple2():
    check_electron_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0.57092, 0.29608, -0.758]), np.array([-0.70841, 0.22864, 0.79589]),
        np.array([0.83984, 0.65053, 0.36087]), np.array([-0.62267, -0.83676, -0.75233]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.19609589]]]]))


def test_electron_repulsion_0_0_0_0_simple3():
    check_electron_repulsion(
        np.array([0.57283]), np.array([1.74713]),
        np.array([0.21032]), np.array([1.60538]),
        np.array([0.82197, 0.73226, -0.98154]), np.array([0.57466, 0.17815, -0.25519]),
        np.array([0.00425, -0.33757, 0.08556]), np.array([-0.38717, 0.66721, 0.40838]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.92553047]]]]))


def test_electron_repulsion_0_0_0_0_simple4():
    check_electron_repulsion(
        np.array([1.35491]), np.array([0.9714]),
        np.array([1.95585]), np.array([1.77853]),
        np.array([0.37263, -0.87382, 0.28078]), np.array([-0.08946, -0.52616, 0.69184]),
        np.array([-0.35128, 0.07017, 0.08193]), np.array([0.14543, -0.29499, -0.09769]),
        np.array([1.61086]),
        np.array([1.19397]),
        np.array([1.8119]),
        np.array([1.55646]),
        0, 0, 0, 0,
        np.array([[[[1.65373353]]]]))


def test_electron_repulsion_0_0_0_0():
    check_electron_repulsion(
        np.array([1.63216, 1.25493, 1.46134, 0.48024]),
        np.array([1.72365, 1.59905, 0.10447, 1.28324]),
        np.array([1.4105, 0.27134, 1.51238, 0.7518]),
        np.array([1.38488, 0.97611, 0.34149, 0.4326]),
        np.array([0.61356, -0.85284, -0.37151]), np.array([-0.63238, -0.81396, 0.40314]),
        np.array([0.29559, 0.60342, 0.18878]), np.array([-0.6893, 0.09175, -0.97283]),
        np.array([0.84965]),
        np.array([1.49169]),
        np.array([1.11046]),
        np.array([0.6665]),
        0, 0, 0, 0,
        np.array([[[[2.97326773]]]]))


def test_electron_repulsion_0_0_0_1():
    check_electron_repulsion(
        np.array([0.74579, 0.93686, 0.39742]), np.array([1.01349, 1.46072, 0.22295]),
        np.array([1.90756, 0.52423, 1.35586]), np.array([0.9655, 0.73539, 0.51017]),
        np.array([0.55177, 0.11232, -0.95152]), np.array([0.79941, 0.80782, 0.02287]),
        np.array([-0.52471, 0.59124, 0.434]), np.array([0.40758, 0.96818, 0.59852]),
        np.array([1.38989]),
        np.array([1.20619]),
        np.array([1.25917]),
        np.array([0.70246, 1.69253, 1.5632]),
        0, 0, 0, 1,
        np.array([[[[0.03331693, -3.27978915, -5.8871596]]]]))


def test_electron_repulsion_0_0_1_0():
    check_electron_repulsion(
        np.array([0.88609, 0.76883, 0.56082]), np.array([1.29216, 0.28671, 1.25389]),
        np.array([1.36987, 0.90792, 0.30511]), np.array([0.57079, 1.98163, 0.66835]),
        np.array([0.7706, 0.99091, -0.21592]), np.array([-0.00566, -0.37522, -0.3936]),
        np.array([0.1527, -0.95347, 0.16682]), np.array([-0.75347, -0.6388, -0.81567]),
        np.array([0.72776]),
        np.array([1.08088]),
        np.array([0.53874, 1.37722, 1.16945]),
        np.array([1.4472]),
        0, 0, 1, 0,
        np.array([[[[0.5110845], [5.86518324], [-1.47878266]]]]))


def test_electron_repulsion_0_0_1_1():
    check_electron_repulsion(
        np.array([0.94138, 0.23708, 1.33464]), np.array([1.89753, 0.54214, 0.80346]),
        np.array([1.04131, 1.6925, 0.81454]), np.array([1.06467, 0.55116, 1.21121]),
        np.array([0.6941, 0.3354, -0.49162]), np.array([0.68756, 0.49975, -0.69756]),
        np.array([0.60432, -0.01449, -0.26057]), np.array([0.35763, -0.04674, -0.78137]),
        np.array([0.75847]),
        np.array([0.57683]),
        np.array([1.61747, 0.59289, 0.93361]),
        np.array([1.38523, 1.77715, 0.8249]),
        0, 0, 1, 1,
        np.array([[[[0.57999607, 0.04732015, 0.00079488],
                    [0.09549513, 0.42707461, 0.03630467],
                    [-0.15902635, -0.25704193, 0.12295133]]]]))


def test_electron_repulsion_0_1_0_0():
    check_electron_repulsion(
        np.array([0.11308, 0.49861, 1.12215]), np.array([0.6186, 1.93501, 1.72751]),
        np.array([0.4644, 0.61371, 1.99408]), np.array([1.98686, 0.49338, 0.88466]),
        np.array([0.31794, 0.18412, 0.89808]), np.array([0.35463, 0.17042, 0.0682]),
        np.array([0.51676, -0.86674, -0.32785]), np.array([-0.03453, -0.05741, -0.86135]),
        np.array([1.84487]),
        np.array([1.17293, 1.02836, 0.50605]),
        np.array([0.54734]),
        np.array([1.55774]),
        0, 1, 0, 0,
        np.array([[[[-2.98984233]], [[-2.16665085]], [[-3.19087757]]]]))


def test_electron_repulsion_0_1_0_1():
    check_electron_repulsion(
        np.array([0.95345, 1.7616, 0.62144]), np.array([0.60537, 0.78954, 0.17662]),
        np.array([1.39946, 1.03161, 1.42837]), np.array([1.05228, 1.80211, 1.37614]),
        np.array([0.18086, -0.0927, -0.36495]), np.array([0.48062, -0.97782, -0.05878]),
        np.array([-0.55927, -0.95238, 0.33122]), np.array([0.17856, 0.06077, 0.62697]),
        np.array([0.9876]),
        np.array([1.39633, 1.30787, 1.80682]),
        np.array([0.93201]),
        np.array([1.21516, 1.84023, 1.59345]),
        0, 1, 0, 1,
        np.array([[[[1.11620596, 0.60061237, 0.36843148]],
                   [[-0.05340867, 0.33119515, -0.70418275]],
                   [[-0.04504112, -1.01394262, 1.17313632]]]]))


def test_electron_repulsion_0_1_1_1():
    check_electron_repulsion(
        np.array([1.60961, 1.48434, 1.09022]), np.array([1.49016, 0.78972, 1.01383]),
        np.array([1.357, 1.6929, 1.46297]), np.array([1.3126, 1.39773, 0.3295]),
        np.array([-0.74441, 0.13168, 0.17287]), np.array([-0.73242, 0.73598, -0.07688]),
        np.array([0.06303, 0.61361, 0.92689]), np.array([0.31395, 0.00081, -0.13425]),
        np.array([1.92653]),
        np.array([0.84324, 1.68215, 0.64055]),
        np.array([1.62317, 1.94784, 1.54325]),
        np.array([0.67873, 0.76053, 0.57816]),
        0, 1, 1, 1,
        np.array([[[[-0.06633908, -0.13761956, -0.03005655],
                    [-0.023407, -0.07813472, -0.03489736],
                    [-0.02263273, -0.20143856, -0.03550443]],
                   [[-0.40044718, -0.35436776, 0.07827812],
                    [-0.39382673, -0.18295174, 0.10845718],
                    [-0.37310311, -0.34400264, 0.05152883]],
                   [[0.07743294, -0.04648822, -0.2043075],
                    [0.03540926, -0.00400861, -0.13446393],
                    [0.02364929, -0.01807209, -0.18079094]]]]))


def test_electron_repulsion_1_0_0_1():
    check_electron_repulsion(
        np.array([0.39834, 1.4798, 1.80662]), np.array([1.9623, 0.88607, 0.93517]),
        np.array([0.46864, 1.1317, 0.67625]), np.array([1.52214, 0.93879, 0.71425]),
        np.array([-0.04796, 0.70504, 0.36481]), np.array([0.40599, 0.97607, 0.64758]),
        np.array([0.66271, -0.64123, -0.17474]), np.array([-0.60087, 0.25093, 0.32664]),
        np.array([0.68301, 1.18047, 1.44482]),
        np.array([0.97181]),
        np.array([1.18315]),
        np.array([0.79184, 1.41932, 1.32812]),
        1, 0, 0, 1,
        np.array([[[[0.16173756, 0.14265052, 0.05405344]]],
                  [[[-0.431925, -0.37295006, -0.1782411]]],
                  [[[-0.17915755, -0.20235955, 0.03526912]]]]))


def test_electron_repulsion_1_1_1_1():
    check_electron_repulsion(
        np.array([0.13992, 0.37329, 0.33259]), np.array([0.64139, 1.73019, 0.13917]),
        np.array([0.44337, 1.28161, 0.3277]), np.array([1.24252, 1.27924, 1.45445]),
        np.array([0.02582, 0.94923, -0.17438]), np.array([-0.81301, 0.086, -0.77236]),
        np.array([-0.67901, 0.6566, -0.45438]), np.array([-0.02669, -0.13942, -0.98892]),
        np.array([1.01729, 0.83942, 1.15976]),
        np.array([1.92943, 1.10829, 0.87557]),
        np.array([0.58667, 0.97031, 1.31261]),
        np.array([1.57111, 0.74218, 0.68171]),
        1, 1, 1, 1,
        np.array([[[[5.38092832, 0.67101024, 0.50643354],
                    [-0.36637823, -0.17128347, 0.00749151],
                    [-0.47015285, -0.00846274, -0.23514519]],
                   [[0.31412053, 1.85552661, -0.05096966],
                    [0.5668773, 0.04019152, -0.05803149],
                    [-0.02195855, 0.00256108, 0.03373068]],
                   [[0.26139911, -0.05908764, 1.34729127],
                    [0.03563575, 0.02599451, 0.0669569],
                    [0.6249628, -0.09012696, -0.02559206]]],
                  [[[-1.30079959, 0.06525516, -0.24130176],
                    [7.90805546, 0.5029288, 1.03164863],
                    [0.22531828, -0.01518479, -0.63472654]],
                   [[0.07758755, -0.30344079, 0.03679751],
                    [0.88274549, 3.43263474, -0.20761467],
                    [0.09249023, 0.10854722, 0.15741632]],
                   [[0.0082139, -0.00382022, -0.24202072],
                    [0.44155444, -0.06437548, 2.40552259],
                    [0.29276089, 0.01725224, 0.05956368]]],
                  [[[-1.45339037, -0.37266055, 0.25844897],
                    [0.41152374, -0.40525461, -0.16607501],
                    [14.23224926, 2.34068558, 0.65653732]],
                   [[-0.00776144, -0.38261119, -0.0073076],
                    [0.28311943, 0.14089539, 0.08426703],
                    [0.91304633, 5.92042353, -0.12886949]],
                   [[0.09807363, 0.06281554, -0.25920407],
                    [0.15636252, 0.10752926, 0.14182457],
                    [1.2142302, -0.38098265, 4.57694241]]]]))


def test_electron_repulsion_0_2_1_0():
    check_electron_repulsion(
        np.array([1.36794, 1.14001, 1.97798]), np.array([1.68538, 0.75019, 0.72741]),
        np.array([1.55248, 0.78842, 1.84644]), np.array([1.73266, 0.46153, 0.63621]),
        np.array([0.05517, 0.27196, -0.98928]), np.array([-0.20526, 0.27314, -0.16208]),
        np.array([-0.00876, -0.47585, 0.88613]), np.array([0.75034, 0.54371, -0.1464]),
        np.array([0.50974]),
        np.array([1.36246, 0.58913, 0.73488, 0.53568, 1.11864, 1.80388]),
        np.array([1.62815, 0.58942, 1.52452]),
        np.array([0.66094]),
        0, 2, 1, 0,
        np.array([[[[0.03940319], [0.05597157], [-0.32990373]],
                   [[-0.00066587], [0.00221213], [-0.00319745]],
                   [[-0.00035194], [-0.00011777], [0.0063613]],
                   [[0.00478058], [0.01592957], [-0.09687372]],
                   [[0.00002574], [-0.00009517], [-0.00166564]],
                   [[0.01578456], [0.05420504], [-0.32175899]]]]))


def test_electron_repulsion_0_2_2_3():
    fn = 'electron_repulsion_0_2_2_3_json'
    result0 = load_json(fn)
    check_electron_repulsion(
        np.array([0.96867, 0.41743, 1.03509]), np.array([1.84594, 0.83035, 1.20242]),
        np.array([0.94861, 0.47292, 0.38655]), np.array([1.3009, 1.10486, 1.4979]),
        np.array([0.10017, 0.21708, 0.08942]), np.array([-0.03049, 0.99486, -0.37959]),
        np.array([-0.7765, 0.53988, 0.25643]), np.array([0.60758, 0.85146, 0.15088]),
        np.array([1.14284]),
        np.array([1.39723, 1.77896, 0.72525, 0.99877, 1.5953, 0.69473]),
        np.array([0.56774, 1.69348, 1.8146, 0.85426, 1.35434, 1.87402]),
        np.array([0.99964, 1.45499, 1.35143, 1.9758, 0.58887, 1.40713, 0.55226, 1.44979,
                  0.57156, 0.71009]),
        0, 2, 2, 3,
        result0)


def test_electron_repulsion_4_3_2_1():
    fn = 'electron_repulsion_4_3_2_1_json'
    result0 = load_json(fn)
    check_electron_repulsion(
        np.array([0.94212, 1.71823, 0.3309]), np.array([0.94854, 0.12816, 0.42016]),
        np.array([0.46046, 0.43321, 1.0587]), np.array([1.0089, 0.52286, 1.83539]),
        np.array([-0.48859, 0.6043, -0.57858]), np.array([0.74567, -0.82555, -0.30631]),
        np.array([-0.5679, -0.08725, 0.7623]), np.array([0.10338, 0.65407, -0.20172]),
        np.array([1.10904, 1.40637, 1.8707, 0.68295, 1.29692, 0.99892, 1.13936, 0.81258,
                  0.50325, 1.27698, 1.81192, 1.43415, 1.1686, 1.38063, 0.61592]),
        np.array([1.19368, 0.75291, 0.63535, 1.22654, 1.32848, 1.17482, 1.74897, 0.93964,
                  1.90303, 1.44528]),
        np.array([1.63343, 1.80498, 1.61313, 0.99992, 1.04505, 1.42297]),
        np.array([1.4825, 1.69421, 1.8635]),
        4, 3, 2, 1,
        result0)


def get_erf_repulsion(alphas0, alphas1, alphas2, alphas3, r0, r1, r2, r3,
                      scales0, scales1, scales2, scales3,
                      shell_type0, shell_type1, shell_type2, shell_type3, mu):
    """Get the short-range damped Erf integrals for a primitive shell.

    Parameters
    ----------
    alpha0, alpha1, alpha2, alpha3 : float
        Exponents of the four primitive shells.
    r0, r1, r2, r3 : np.ndarray, shape=(3,), dtype=float
        Cartesian coordinates of the centers of the four primitive shells.
    scales0, scales1, scales2, scales3 : np.ndarray, dtype=float
        Normalization prefactors for the Gaussian shells.
    shell_type0, shell_type1, shell_type2, shell_type3 : int
        Shell types of the four primitive shells.
    mu : float
        The range-separation parameters.
    """
    max_shell_type = 4
    gb4i = _GB4ErfIntegralLibInt(max_shell_type, mu)

    nbasis0 = _get_shell_nbasis(shell_type0)
    nbasis1 = _get_shell_nbasis(shell_type1)
    nbasis2 = _get_shell_nbasis(shell_type2)
    nbasis3 = _get_shell_nbasis(shell_type3)
    # Clear the working memory
    gb4i.reset(shell_type0, shell_type1, shell_type2, shell_type3, r0, r1, r2, r3)
    # Add a few contributions:
    for alpha0, alpha1, alpha2, alpha3 in zip(alphas0, alphas1, alphas2, alphas3):
        gb4i.add(1.0, alpha0, alpha1, alpha2, alpha3, scales0, scales1, scales2, scales3)
    return gb4i.get_work(nbasis0, nbasis1, nbasis2, nbasis3)


def check_erf_repulsion(alphas0, alphas1, alphas2, alphas3, r0, r1, r2, r3, scales0,
                        scales1, scales2, scales3, shell_type0, shell_type1, shell_type2,
                        shell_type3, result0, mu):
    """Compare output from HORTON Erf integrals with reference data.

    The reference data was generated with a Mathematica script of Julien Toulouse and
    Andreas Savin.

    Parameters
    ----------
    alpha0, alpha1, alpha2, alpha3 : float
        Exponents of the four primitive shells.
    r0, r1, r2, r3 : np.ndarray, shape=(3,), dtype=float
        Cartesian coordinates of the centers of the four primitive shells.
    scales0, scales1, scales2, scales3 : np.ndarray, dtype=float
        Normalization prefactors for the Gaussian shells.
    shell_type0, shell_type1, shell_type2, shell_type3 : int
        Shell types of the four primitive shells.
    result0 : np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis), dtype=float
        The expected result.
    mu : float
        The range-separation parameters.
    """
    result1 = get_erf_repulsion(alphas0, alphas1, alphas2, alphas3, r0, r1, r2, r3,
                                scales0, scales1, scales2, scales3, shell_type0,
                                shell_type1, shell_type2, shell_type3, mu)
    assert abs(result1 - result0).max() < 3e-7


def test_erf_repulsion_0_0_0_0_simple0():
    check_erf_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[1.25667419]]]]), 0.3)


def test_erf_repulsion_0_0_0_0_simple1():
    check_erf_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[1.16018914]]]]), 0.3)


def test_erf_repulsion_0_0_0_0_simple2():
    check_erf_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0.57092, 0.29608, -0.758]), np.array([-0.70841, 0.22864, 0.79589]),
        np.array([0.83984, 0.65053, 0.36087]), np.array([-0.62267, -0.83676, -0.75233]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.09691428]]]]), 0.3)


def test_erf_repulsion_0_0_0_0_simple3():
    check_erf_repulsion(
        np.array([0.57283]), np.array([1.74713]),
        np.array([0.21032]), np.array([1.60538]),
        np.array([0.82197, 0.73226, -0.98154]), np.array([0.57466, 0.17815, -0.25519]),
        np.array([0.00425, -0.33757, 0.08556]), np.array([-0.38717, 0.66721, 0.40838]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.807980035]]]]), 1.2)


def test_erf_repulsion_0_0_0_0_simple4():
    check_erf_repulsion(
        np.array([1.35491]), np.array([0.9714]),
        np.array([1.95585]), np.array([1.77853]),
        np.array([0.37263, -0.87382, 0.28078]), np.array([-0.08946, -0.52616, 0.69184]),
        np.array([-0.35128, 0.07017, 0.08193]), np.array([0.14543, -0.29499, -0.09769]),
        np.array([1.61086]),
        np.array([1.19397]),
        np.array([1.8119]),
        np.array([1.55646]),
        0, 0, 0, 0,
        np.array([[[[1.16218348]]]]), 1.2)


def test_erf_repulsion_0_0_0_0():
    check_erf_repulsion(
        np.array([1.63216, 1.25493, 1.46134, 0.48024]),
        np.array([1.72365, 1.59905, 0.10447, 1.28324]),
        np.array([1.4105, 0.27134, 1.51238, 0.7518]),
        np.array([1.38488, 0.97611, 0.34149, 0.4326]),
        np.array([0.61356, -0.85284, -0.37151]), np.array([-0.63238, -0.81396, 0.40314]),
        np.array([0.29559, 0.60342, 0.18878]), np.array([-0.6893, 0.09175, -0.97283]),
        np.array([0.84965]),
        np.array([1.49169]),
        np.array([1.11046]),
        np.array([0.6665]),
        0, 0, 0, 0,
        np.array([[[[2.97326773]]]]), 10000)


def test_erf_repulsion_0_0_0_1():
    check_erf_repulsion(
        np.array([0.74579, 0.93686, 0.39742]), np.array([1.01349, 1.46072, 0.22295]),
        np.array([1.90756, 0.52423, 1.35586]), np.array([0.9655, 0.73539, 0.51017]),
        np.array([0.55177, 0.11232, -0.95152]), np.array([0.79941, 0.80782, 0.02287]),
        np.array([-0.52471, 0.59124, 0.434]), np.array([0.40758, 0.96818, 0.59852]),
        np.array([1.38989]),
        np.array([1.20619]),
        np.array([1.25917]),
        np.array([0.70246, 1.69253, 1.5632]),
        0, 0, 0, 1,
        np.array([[[[0.03331693, -3.27978915, -5.8871596]]]]), 10000)


def test_erf_repulsion_0_0_1_0():
    check_erf_repulsion(
        np.array([0.88609, 0.76883, 0.56082]), np.array([1.29216, 0.28671, 1.25389]),
        np.array([1.36987, 0.90792, 0.30511]), np.array([0.57079, 1.98163, 0.66835]),
        np.array([0.7706, 0.99091, -0.21592]), np.array([-0.00566, -0.37522, -0.3936]),
        np.array([0.1527, -0.95347, 0.16682]), np.array([-0.75347, -0.6388, -0.81567]),
        np.array([0.72776]),
        np.array([1.08088]),
        np.array([0.53874, 1.37722, 1.16945]),
        np.array([1.4472]),
        0, 0, 1, 0,
        np.array([[[[0.5110845], [5.86518324], [-1.47878266]]]]), 10000)


def test_erf_repulsion_0_0_1_1():
    check_erf_repulsion(
        np.array([0.94138, 0.23708, 1.33464]), np.array([1.89753, 0.54214, 0.80346]),
        np.array([1.04131, 1.6925, 0.81454]), np.array([1.06467, 0.55116, 1.21121]),
        np.array([0.6941, 0.3354, -0.49162]), np.array([0.68756, 0.49975, -0.69756]),
        np.array([0.60432, -0.01449, -0.26057]), np.array([0.35763, -0.04674, -0.78137]),
        np.array([0.75847]),
        np.array([0.57683]),
        np.array([1.61747, 0.59289, 0.93361]),
        np.array([1.38523, 1.77715, 0.8249]),
        0, 0, 1, 1,
        np.array([[[[0.57999607, 0.04732015, 0.00079488],
                    [0.09549513, 0.42707461, 0.03630467],
                    [-0.15902635, -0.25704193, 0.12295133]]]]),
        10000)


def test_erf_repulsion_0_1_0_0():
    check_erf_repulsion(
        np.array([0.11308, 0.49861, 1.12215]), np.array([0.6186, 1.93501, 1.72751]),
        np.array([0.4644, 0.61371, 1.99408]), np.array([1.98686, 0.49338, 0.88466]),
        np.array([0.31794, 0.18412, 0.89808]), np.array([0.35463, 0.17042, 0.0682]),
        np.array([0.51676, -0.86674, -0.32785]), np.array([-0.03453, -0.05741, -0.86135]),
        np.array([1.84487]),
        np.array([1.17293, 1.02836, 0.50605]),
        np.array([0.54734]),
        np.array([1.55774]),
        0, 1, 0, 0,
        np.array([[[[-2.98984233]], [[-2.16665085]], [[-3.19087757]]]]), 10000)


def test_erf_repulsion_0_1_0_1():
    check_erf_repulsion(
        np.array([0.95345, 1.7616, 0.62144]), np.array([0.60537, 0.78954, 0.17662]),
        np.array([1.39946, 1.03161, 1.42837]), np.array([1.05228, 1.80211, 1.37614]),
        np.array([0.18086, -0.0927, -0.36495]), np.array([0.48062, -0.97782, -0.05878]),
        np.array([-0.55927, -0.95238, 0.33122]), np.array([0.17856, 0.06077, 0.62697]),
        np.array([0.9876]),
        np.array([1.39633, 1.30787, 1.80682]),
        np.array([0.93201]),
        np.array([1.21516, 1.84023, 1.59345]),
        0, 1, 0, 1,
        np.array([[[[1.11620596, 0.60061237, 0.36843148]],
                   [[-0.05340867, 0.33119515, -0.70418275]],
                   [[-0.04504112, -1.01394262, 1.17313632]]]]),
        10000)


def test_erf_repulsion_0_1_1_1():
    check_erf_repulsion(
        np.array([1.60961, 1.48434, 1.09022]), np.array([1.49016, 0.78972, 1.01383]),
        np.array([1.357, 1.6929, 1.46297]), np.array([1.3126, 1.39773, 0.3295]),
        np.array([-0.74441, 0.13168, 0.17287]), np.array([-0.73242, 0.73598, -0.07688]),
        np.array([0.06303, 0.61361, 0.92689]), np.array([0.31395, 0.00081, -0.13425]),
        np.array([1.92653]),
        np.array([0.84324, 1.68215, 0.64055]),
        np.array([1.62317, 1.94784, 1.54325]),
        np.array([0.67873, 0.76053, 0.57816]),
        0, 1, 1, 1,
        np.array([[[[-0.06633908, -0.13761956, -0.03005655],
                    [-0.023407, -0.07813472, -0.03489736],
                    [-0.02263273, -0.20143856, -0.03550443]],
                   [[-0.40044718, -0.35436776, 0.07827812],
                    [-0.39382673, -0.18295174, 0.10845718],
                    [-0.37310311, -0.34400264, 0.05152883]],
                   [[0.07743294, -0.04648822, -0.2043075],
                    [0.03540926, -0.00400861, -0.13446393],
                    [0.02364929, -0.01807209, -0.18079094]]]]),
        10000)


def test_erf_repulsion_1_0_0_1():
    check_erf_repulsion(
        np.array([0.39834, 1.4798, 1.80662]), np.array([1.9623, 0.88607, 0.93517]),
        np.array([0.46864, 1.1317, 0.67625]), np.array([1.52214, 0.93879, 0.71425]),
        np.array([-0.04796, 0.70504, 0.36481]), np.array([0.40599, 0.97607, 0.64758]),
        np.array([0.66271, -0.64123, -0.17474]), np.array([-0.60087, 0.25093, 0.32664]),
        np.array([0.68301, 1.18047, 1.44482]),
        np.array([0.97181]),
        np.array([1.18315]),
        np.array([0.79184, 1.41932, 1.32812]),
        1, 0, 0, 1,
        np.array([[[[0.16173756, 0.14265052, 0.05405344]]],
                  [[[-0.431925, -0.37295006, -0.1782411]]],
                  [[[-0.17915755, -0.20235955, 0.03526912]]]]),
        10000)


def test_erf_repulsion_1_1_1_1():
    check_erf_repulsion(
        np.array([0.13992, 0.37329, 0.33259]), np.array([0.64139, 1.73019, 0.13917]),
        np.array([0.44337, 1.28161, 0.3277]), np.array([1.24252, 1.27924, 1.45445]),
        np.array([0.02582, 0.94923, -0.17438]), np.array([-0.81301, 0.086, -0.77236]),
        np.array([-0.67901, 0.6566, -0.45438]), np.array([-0.02669, -0.13942, -0.98892]),
        np.array([1.01729, 0.83942, 1.15976]),
        np.array([1.92943, 1.10829, 0.87557]),
        np.array([0.58667, 0.97031, 1.31261]),
        np.array([1.57111, 0.74218, 0.68171]),
        1, 1, 1, 1,
        np.array([[[[5.38092832, 0.67101024, 0.50643354],
                    [-0.36637823, -0.17128347, 0.00749151],
                    [-0.47015285, -0.00846274, -0.23514519]],
                   [[0.31412053, 1.85552661, -0.05096966],
                    [0.5668773, 0.04019152, -0.05803149],
                    [-0.02195855, 0.00256108, 0.03373068]],
                   [[0.26139911, -0.05908764, 1.34729127],
                    [0.03563575, 0.02599451, 0.0669569],
                    [0.6249628, -0.09012696, -0.02559206]]],
                  [[[-1.30079959, 0.06525516, -0.24130176],
                    [7.90805546, 0.5029288, 1.03164863],
                    [0.22531828, -0.01518479, -0.63472654]],
                   [[0.07758755, -0.30344079, 0.03679751],
                    [0.88274549, 3.43263474, -0.20761467],
                    [0.09249023, 0.10854722, 0.15741632]],
                   [[0.0082139, -0.00382022, -0.24202072],
                    [0.44155444, -0.06437548, 2.40552259],
                    [0.29276089, 0.01725224, 0.05956368]]],
                  [[[-1.45339037, -0.37266055, 0.25844897],
                    [0.41152374, -0.40525461, -0.16607501],
                    [14.23224926, 2.34068558, 0.65653732]],
                   [[-0.00776144, -0.38261119, -0.0073076],
                    [0.28311943, 0.14089539, 0.08426703],
                    [0.91304633, 5.92042353, -0.12886949]],
                   [[0.09807363, 0.06281554, -0.25920407],
                    [0.15636252, 0.10752926, 0.14182457],
                    [1.2142302, -0.38098265, 4.57694241]]]]),
        10000)


def test_erf_repulsion_0_2_1_0():
    check_erf_repulsion(
        np.array([1.36794, 1.14001, 1.97798]), np.array([1.68538, 0.75019, 0.72741]),
        np.array([1.55248, 0.78842, 1.84644]), np.array([1.73266, 0.46153, 0.63621]),
        np.array([0.05517, 0.27196, -0.98928]), np.array([-0.20526, 0.27314, -0.16208]),
        np.array([-0.00876, -0.47585, 0.88613]), np.array([0.75034, 0.54371, -0.1464]),
        np.array([0.50974]),
        np.array([1.36246, 0.58913, 0.73488, 0.53568, 1.11864, 1.80388]),
        np.array([1.62815, 0.58942, 1.52452]),
        np.array([0.66094]),
        0, 2, 1, 0,
        np.array([[[[0.03940319], [0.05597157], [-0.32990373]],
                   [[-0.00066587], [0.00221213], [-0.00319745]],
                   [[-0.00035194], [-0.00011777], [0.0063613]],
                   [[0.00478058], [0.01592957], [-0.09687372]],
                   [[0.00002574], [-0.00009517], [-0.00166564]],
                   [[0.01578456], [0.05420504], [-0.32175899]]]]),
        10000)


def test_erf_repulsion_0_2_2_3():
    fn = 'erf_repulsion_0_2_2_3_json'
    result0 = load_json(fn)
    check_erf_repulsion(
        np.array([0.96867, 0.41743, 1.03509]), np.array([1.84594, 0.83035, 1.20242]),
        np.array([0.94861, 0.47292, 0.38655]), np.array([1.3009, 1.10486, 1.4979]),
        np.array([0.10017, 0.21708, 0.08942]), np.array([-0.03049, 0.99486, -0.37959]),
        np.array([-0.7765, 0.53988, 0.25643]), np.array([0.60758, 0.85146, 0.15088]),
        np.array([1.14284]),
        np.array([1.39723, 1.77896, 0.72525, 0.99877, 1.5953, 0.69473]),
        np.array([0.56774, 1.69348, 1.8146, 0.85426, 1.35434, 1.87402]),
        np.array([0.99964, 1.45499, 1.35143, 1.9758, 0.58887, 1.40713, 0.55226,
                  1.44979, 0.57156, 0.71009]),
        0, 2, 2, 3,
        result0, 10000)


def test_erf_repulsion_4_3_2_1():
    fn = 'erf_repulsion_4_3_2_1_json'
    result0 = load_json(fn)
    check_erf_repulsion(
        np.array([0.94212, 1.71823, 0.3309]), np.array([0.94854, 0.12816, 0.42016]),
        np.array([0.46046, 0.43321, 1.0587]), np.array([1.0089, 0.52286, 1.83539]),
        np.array([-0.48859, 0.6043, -0.57858]), np.array([0.74567, -0.82555, -0.30631]),
        np.array([-0.5679, -0.08725, 0.7623]), np.array([0.10338, 0.65407, -0.20172]),
        np.array([1.10904, 1.40637, 1.8707, 0.68295, 1.29692, 0.99892, 1.13936, 0.81258,
                  0.50325, 1.27698, 1.81192, 1.43415, 1.1686, 1.38063, 0.61592]),
        np.array([1.19368, 0.75291, 0.63535, 1.22654, 1.32848, 1.17482, 1.74897, 0.93964,
                  1.90303, 1.44528]),
        np.array([1.63343, 1.80498, 1.61313, 0.99992, 1.04505, 1.42297]),
        np.array([1.4825, 1.69421, 1.8635]),
        4, 3, 2, 1,
        result0, 10000)


# TODO: Move to glue code test?
# def test_erf_repulsion_h2_sto3g():
#     mol2 = IOData.from_file('FCIDUMP.molpro.h2-erf'))
#
#     mol = IOData(title='h2')
#     mol['coordinates'] = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])
#     mol['numbers'] = np.array([1, 1])
#
#     # Create a Gaussian basis set
#     obasis = get_gobasis(mol['coordinates'], mol['numbers'], 'sto-3g')
#
#     # Compute Gaussian integrals
#     olp = obasis.compute_overlap()
#     kin = obasis.compute_kinetic()
#     na = obasis.compute_nuclear_attraction(mol['coordinates'], mol['pseudo_numbers'])
#     er = obasis.compute_erf_repulsion(2.25)
#
#     # Create alpha orbitals
#     orb_alpha = Orbitals(obasis.nbasis)
#
#     # Initial guess
#     core = kin + na
#     guess_core_hamiltonian(olp, core, orb_alpha)
#     mol.orb_alpha = orb_alpha
#
#     # Transform orbitals
#     two_mo = transform_integrals(core, er, 'tensordot', mol.orb_alpha)[1][0]
#     assert abs(mol2.two_mo - two_mo).max() < 1e-10


def check_gauss_repulsion(alphas0, alphas1, alphas2, alphas3, r0, r1, r2, r3, scales0,
                          scales1, scales2, scales3, shell_type0, shell_type1,
                          shell_type2, shell_type3, result0, c, alpha):
    """Compare output from HORTON Gauss 4-center integrals with reference data.

    The reference data was generated with a Mathematica script of Julien Toulouse and
    Andreas Savin.

    Parameters
    ----------
    alpha0, alpha1, alpha2, alpha3 : float
        Exponents of the four primitive shells.
    r0, r1, r2, r3 : np.ndarray, shape=(3,), dtype=float
        Cartesian coordinates of the centers of the four primitive shells.
    scales0, scales1, scales2, scales3 : np.ndarray, dtype=float
        Normalization prefactors for the Gaussian shells.
    shell_type0, shell_type1, shell_type2, shell_type3 : int
        Shell types of the four primitive shells.
    result0 : np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis), dtype=float
        The expected result.
    c : float
        Coefficient of the gaussian.
    alpha : float
        Exponential parameter of the gaussian.
    """
    max_shell_type = 4
    max_nbasis = _get_shell_nbasis(max_shell_type)
    gb4i = _GB4GaussIntegralLibInt(max_shell_type, c, alpha)
    assert gb4i.max_nbasis == max_nbasis
    assert gb4i.alpha == alpha
    assert gb4i.c == c
    assert gb4i.nwork == max_nbasis ** 4

    nbasis0 = _get_shell_nbasis(shell_type0)
    nbasis1 = _get_shell_nbasis(shell_type1)
    nbasis2 = _get_shell_nbasis(shell_type2)
    nbasis3 = _get_shell_nbasis(shell_type3)
    assert result0.shape == (nbasis0, nbasis1, nbasis2, nbasis3)
    # Clear the working memory
    gb4i.reset(shell_type0, shell_type1, shell_type2, shell_type3, r0, r1, r2, r3)
    # Add a few contributions:
    for alpha0, alpha1, alpha2, alpha3 in zip(alphas0, alphas1, alphas2, alphas3):
        gb4i.add(1.0, alpha0, alpha1, alpha2, alpha3, scales0, scales1, scales2, scales3)
    result1 = gb4i.get_work(nbasis0, nbasis1, nbasis2, nbasis3)
    assert abs(result1 - result0).max() < 3e-7


def test_gauss_repulsion_0_0_0_0_simple0():
    check_gauss_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[-2.45402796]]]]), -2.256758334191, 4. / 3.)


def test_gauss_repulsion_0_0_0_0_simple1():
    check_gauss_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[1.1063809195]]]]), -0.5 * -2.25675833419102, 1.5)


def test_gauss_repulsion_0_0_0_0_simple2():
    check_gauss_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[-0.44195157099]]]]), -2.25675833419102, 4. / 3.)


def test_gauss_repulsion_0_0_0_0_simple3():
    check_gauss_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0.57092, 0.29608, -0.758]), np.array([-0.70841, 0.22864, 0.79589]),
        np.array([0.83984, 0.65053, 0.36087]), np.array([-0.62267, -0.83676, -0.75233]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[-0.0476487494]]]]), -2.25675833419102, 4. / 3.)


def test_gauss_repulsion_0_0_0_0_simple4():
    check_gauss_repulsion(
        np.array([0.57283]), np.array([1.74713]),
        np.array([0.21032]), np.array([1.60538]),
        np.array([0.82197, 0.73226, -0.98154]), np.array([0.57466, 0.17815, -0.25519]),
        np.array([0.00425, -0.33757, 0.08556]), np.array([-0.38717, 0.66721, 0.40838]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[-0.352255229]]]]), -2.25675833419102, 4. / 3.)


def test_gauss_repulsion_0_0_0_0_simple5():
    check_gauss_repulsion(
        np.array([1.35491]), np.array([0.9714]),
        np.array([1.95585]), np.array([1.77853]),
        np.array([0.37263, -0.87382, 0.28078]), np.array([-0.08946, -0.52616, 0.69184]),
        np.array([-0.35128, 0.07017, 0.08193]), np.array([0.14543, -0.29499, -0.09769]),
        np.array([1.61086]),
        np.array([1.19397]),
        np.array([1.8119]),
        np.array([1.55646]),
        0, 0, 0, 0,
        np.array([[[[-1.03673858]]]]), -2.25675833419102, 4. / 3.)


def test_gauss_repulsion_0_0_0_1():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([0.74579, 0.93686, 0.39742]), np.array([1.01349, 1.46072, 0.22295]),
        np.array([1.90756, 0.52423, 1.35586]), np.array([0.9655, 0.73539, 0.51017]),
        np.array([0.55177, 0.11232, -0.95152]), np.array([0.79941, 0.80782, 0.02287]),
        np.array([-0.52471, 0.59124, 0.434]), np.array([0.40758, 0.96818, 0.59852]),
        np.array([1.38989]),
        np.array([1.20619]),
        np.array([1.25917]),
        np.array([0.70246, 1.69253, 1.5632]),
        0, 0, 0, 1, 2.0)
    e2 = get_erf_repulsion(
        np.array([0.74579, 0.93686, 0.39742]), np.array([1.01349, 1.46072, 0.22295]),
        np.array([1.90756, 0.52423, 1.35586]), np.array([0.9655, 0.73539, 0.51017]),
        np.array([0.55177, 0.11232, -0.95152]), np.array([0.79941, 0.80782, 0.02287]),
        np.array([-0.52471, 0.59124, 0.434]), np.array([0.40758, 0.96818, 0.59852]),
        np.array([1.38989]),
        np.array([1.20619]),
        np.array([1.25917]),
        np.array([0.70246, 1.69253, 1.5632]),
        0, 0, 0, 1, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([0.74579, 0.93686, 0.39742]), np.array([1.01349, 1.46072, 0.22295]),
        np.array([1.90756, 0.52423, 1.35586]), np.array([0.9655, 0.73539, 0.51017]),
        np.array([0.55177, 0.11232, -0.95152]), np.array([0.79941, 0.80782, 0.02287]),
        np.array([-0.52471, 0.59124, 0.434]), np.array([0.40758, 0.96818, 0.59852]),
        np.array([1.38989]),
        np.array([1.20619]),
        np.array([1.25917]),
        np.array([0.70246, 1.69253, 1.5632]),
        0, 0, 0, 1, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_0_0_1_0():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([0.88609, 0.76883, 0.56082]), np.array([1.29216, 0.28671, 1.25389]),
        np.array([1.36987, 0.90792, 0.30511]), np.array([0.57079, 1.98163, 0.66835]),
        np.array([0.7706, 0.99091, -0.21592]), np.array([-0.00566, -0.37522, -0.3936]),
        np.array([0.1527, -0.95347, 0.16682]), np.array([-0.75347, -0.6388, -0.81567]),
        np.array([0.72776]),
        np.array([1.08088]),
        np.array([0.53874, 1.37722, 1.16945]),
        np.array([1.4472]),
        0, 0, 1, 0, 2.0)
    e2 = get_erf_repulsion(
        np.array([0.88609, 0.76883, 0.56082]), np.array([1.29216, 0.28671, 1.25389]),
        np.array([1.36987, 0.90792, 0.30511]), np.array([0.57079, 1.98163, 0.66835]),
        np.array([0.7706, 0.99091, -0.21592]), np.array([-0.00566, -0.37522, -0.3936]),
        np.array([0.1527, -0.95347, 0.16682]), np.array([-0.75347, -0.6388, -0.81567]),
        np.array([0.72776]),
        np.array([1.08088]),
        np.array([0.53874, 1.37722, 1.16945]),
        np.array([1.4472]),
        0, 0, 1, 0, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([0.88609, 0.76883, 0.56082]), np.array([1.29216, 0.28671, 1.25389]),
        np.array([1.36987, 0.90792, 0.30511]), np.array([0.57079, 1.98163, 0.66835]),
        np.array([0.7706, 0.99091, -0.21592]), np.array([-0.00566, -0.37522, -0.3936]),
        np.array([0.1527, -0.95347, 0.16682]), np.array([-0.75347, -0.6388, -0.81567]),
        np.array([0.72776]),
        np.array([1.08088]),
        np.array([0.53874, 1.37722, 1.16945]),
        np.array([1.4472]),
        0, 0, 1, 0, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_0_0_1_1():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([0.94138, 0.23708, 1.33464]), np.array([1.89753, 0.54214, 0.80346]),
        np.array([1.04131, 1.6925, 0.81454]), np.array([1.06467, 0.55116, 1.21121]),
        np.array([0.6941, 0.3354, -0.49162]), np.array([0.68756, 0.49975, -0.69756]),
        np.array([0.60432, -0.01449, -0.26057]), np.array([0.35763, -0.04674, -0.78137]),
        np.array([0.75847]),
        np.array([0.57683]),
        np.array([1.61747, 0.59289, 0.93361]),
        np.array([1.38523, 1.77715, 0.8249]),
        0, 0, 1, 1, 2.0)
    e2 = get_erf_repulsion(
        np.array([0.94138, 0.23708, 1.33464]), np.array([1.89753, 0.54214, 0.80346]),
        np.array([1.04131, 1.6925, 0.81454]), np.array([1.06467, 0.55116, 1.21121]),
        np.array([0.6941, 0.3354, -0.49162]), np.array([0.68756, 0.49975, -0.69756]),
        np.array([0.60432, -0.01449, -0.26057]), np.array([0.35763, -0.04674, -0.78137]),
        np.array([0.75847]),
        np.array([0.57683]),
        np.array([1.61747, 0.59289, 0.93361]),
        np.array([1.38523, 1.77715, 0.8249]),
        0, 0, 1, 1, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([0.94138, 0.23708, 1.33464]), np.array([1.89753, 0.54214, 0.80346]),
        np.array([1.04131, 1.6925, 0.81454]), np.array([1.06467, 0.55116, 1.21121]),
        np.array([0.6941, 0.3354, -0.49162]), np.array([0.68756, 0.49975, -0.69756]),
        np.array([0.60432, -0.01449, -0.26057]), np.array([0.35763, -0.04674, -0.78137]),
        np.array([0.75847]),
        np.array([0.57683]),
        np.array([1.61747, 0.59289, 0.93361]),
        np.array([1.38523, 1.77715, 0.8249]),
        0, 0, 1, 1, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_0_1_0_0():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([0.11308, 0.49861, 1.12215]), np.array([0.6186, 1.93501, 1.72751]),
        np.array([0.4644, 0.61371, 1.99408]), np.array([1.98686, 0.49338, 0.88466]),
        np.array([0.31794, 0.18412, 0.89808]), np.array([0.35463, 0.17042, 0.0682]),
        np.array([0.51676, -0.86674, -0.32785]), np.array([-0.03453, -0.05741, -0.86135]),
        np.array([1.84487]),
        np.array([1.17293, 1.02836, 0.50605]),
        np.array([0.54734]),
        np.array([1.55774]),
        0, 1, 0, 0, 2.0)
    e2 = get_erf_repulsion(
        np.array([0.11308, 0.49861, 1.12215]), np.array([0.6186, 1.93501, 1.72751]),
        np.array([0.4644, 0.61371, 1.99408]), np.array([1.98686, 0.49338, 0.88466]),
        np.array([0.31794, 0.18412, 0.89808]), np.array([0.35463, 0.17042, 0.0682]),
        np.array([0.51676, -0.86674, -0.32785]), np.array([-0.03453, -0.05741, -0.86135]),
        np.array([1.84487]),
        np.array([1.17293, 1.02836, 0.50605]),
        np.array([0.54734]),
        np.array([1.55774]),
        0, 1, 0, 0, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([0.11308, 0.49861, 1.12215]), np.array([0.6186, 1.93501, 1.72751]),
        np.array([0.4644, 0.61371, 1.99408]), np.array([1.98686, 0.49338, 0.88466]),
        np.array([0.31794, 0.18412, 0.89808]), np.array([0.35463, 0.17042, 0.0682]),
        np.array([0.51676, -0.86674, -0.32785]), np.array([-0.03453, -0.05741, -0.86135]),
        np.array([1.84487]),
        np.array([1.17293, 1.02836, 0.50605]),
        np.array([0.54734]),
        np.array([1.55774]),
        0, 1, 0, 0, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_0_1_0_1():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([0.95345, 1.7616, 0.62144]), np.array([0.60537, 0.78954, 0.17662]),
        np.array([1.39946, 1.03161, 1.42837]), np.array([1.05228, 1.80211, 1.37614]),
        np.array([0.18086, -0.0927, -0.36495]), np.array([0.48062, -0.97782, -0.05878]),
        np.array([-0.55927, -0.95238, 0.33122]), np.array([0.17856, 0.06077, 0.62697]),
        np.array([0.9876]),
        np.array([1.39633, 1.30787, 1.80682]),
        np.array([0.93201]),
        np.array([1.21516, 1.84023, 1.59345]),
        0, 1, 0, 1, 2.0)
    e2 = get_erf_repulsion(
        np.array([0.95345, 1.7616, 0.62144]), np.array([0.60537, 0.78954, 0.17662]),
        np.array([1.39946, 1.03161, 1.42837]), np.array([1.05228, 1.80211, 1.37614]),
        np.array([0.18086, -0.0927, -0.36495]), np.array([0.48062, -0.97782, -0.05878]),
        np.array([-0.55927, -0.95238, 0.33122]), np.array([0.17856, 0.06077, 0.62697]),
        np.array([0.9876]),
        np.array([1.39633, 1.30787, 1.80682]),
        np.array([0.93201]),
        np.array([1.21516, 1.84023, 1.59345]),
        0, 1, 0, 1, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([0.95345, 1.7616, 0.62144]), np.array([0.60537, 0.78954, 0.17662]),
        np.array([1.39946, 1.03161, 1.42837]), np.array([1.05228, 1.80211, 1.37614]),
        np.array([0.18086, -0.0927, -0.36495]), np.array([0.48062, -0.97782, -0.05878]),
        np.array([-0.55927, -0.95238, 0.33122]), np.array([0.17856, 0.06077, 0.62697]),
        np.array([0.9876]),
        np.array([1.39633, 1.30787, 1.80682]),
        np.array([0.93201]),
        np.array([1.21516, 1.84023, 1.59345]),
        0, 1, 0, 1, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_0_1_1_1():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([1.60961, 1.48434, 1.09022]), np.array([1.49016, 0.78972, 1.01383]),
        np.array([1.357, 1.6929, 1.46297]), np.array([1.3126, 1.39773, 0.3295]),
        np.array([-0.74441, 0.13168, 0.17287]), np.array([-0.73242, 0.73598, -0.07688]),
        np.array([0.06303, 0.61361, 0.92689]), np.array([0.31395, 0.00081, -0.13425]),
        np.array([1.92653]),
        np.array([0.84324, 1.68215, 0.64055]),
        np.array([1.62317, 1.94784, 1.54325]),
        np.array([0.67873, 0.76053, 0.57816]),
        0, 1, 1, 1, 2.0)
    e2 = get_erf_repulsion(
        np.array([1.60961, 1.48434, 1.09022]), np.array([1.49016, 0.78972, 1.01383]),
        np.array([1.357, 1.6929, 1.46297]), np.array([1.3126, 1.39773, 0.3295]),
        np.array([-0.74441, 0.13168, 0.17287]), np.array([-0.73242, 0.73598, -0.07688]),
        np.array([0.06303, 0.61361, 0.92689]), np.array([0.31395, 0.00081, -0.13425]),
        np.array([1.92653]),
        np.array([0.84324, 1.68215, 0.64055]),
        np.array([1.62317, 1.94784, 1.54325]),
        np.array([0.67873, 0.76053, 0.57816]),
        0, 1, 1, 1, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([1.60961, 1.48434, 1.09022]), np.array([1.49016, 0.78972, 1.01383]),
        np.array([1.357, 1.6929, 1.46297]), np.array([1.3126, 1.39773, 0.3295]),
        np.array([-0.74441, 0.13168, 0.17287]), np.array([-0.73242, 0.73598, -0.07688]),
        np.array([0.06303, 0.61361, 0.92689]), np.array([0.31395, 0.00081, -0.13425]),
        np.array([1.92653]),
        np.array([0.84324, 1.68215, 0.64055]),
        np.array([1.62317, 1.94784, 1.54325]),
        np.array([0.67873, 0.76053, 0.57816]),
        0, 1, 1, 1, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_1_0_0_1():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([0.39834, 1.4798, 1.80662]), np.array([1.9623, 0.88607, 0.93517]),
        np.array([0.46864, 1.1317, 0.67625]), np.array([1.52214, 0.93879, 0.71425]),
        np.array([-0.04796, 0.70504, 0.36481]), np.array([0.40599, 0.97607, 0.64758]),
        np.array([0.66271, -0.64123, -0.17474]), np.array([-0.60087, 0.25093, 0.32664]),
        np.array([0.68301, 1.18047, 1.44482]),
        np.array([0.97181]),
        np.array([1.18315]),
        np.array([0.79184, 1.41932, 1.32812]),
        1, 0, 0, 1, 2.0)
    e2 = get_erf_repulsion(
        np.array([0.39834, 1.4798, 1.80662]), np.array([1.9623, 0.88607, 0.93517]),
        np.array([0.46864, 1.1317, 0.67625]), np.array([1.52214, 0.93879, 0.71425]),
        np.array([-0.04796, 0.70504, 0.36481]), np.array([0.40599, 0.97607, 0.64758]),
        np.array([0.66271, -0.64123, -0.17474]), np.array([-0.60087, 0.25093, 0.32664]),
        np.array([0.68301, 1.18047, 1.44482]),
        np.array([0.97181]),
        np.array([1.18315]),
        np.array([0.79184, 1.41932, 1.32812]),
        1, 0, 0, 1, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([0.39834, 1.4798, 1.80662]), np.array([1.9623, 0.88607, 0.93517]),
        np.array([0.46864, 1.1317, 0.67625]), np.array([1.52214, 0.93879, 0.71425]),
        np.array([-0.04796, 0.70504, 0.36481]), np.array([0.40599, 0.97607, 0.64758]),
        np.array([0.66271, -0.64123, -0.17474]), np.array([-0.60087, 0.25093, 0.32664]),
        np.array([0.68301, 1.18047, 1.44482]),
        np.array([0.97181]),
        np.array([1.18315]),
        np.array([0.79184, 1.41932, 1.32812]),
        1, 0, 0, 1, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_1_1_1_1():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([0.13992, 0.37329, 0.33259]), np.array([0.64139, 1.73019, 0.13917]),
        np.array([0.44337, 1.28161, 0.3277]), np.array([1.24252, 1.27924, 1.45445]),
        np.array([0.02582, 0.94923, -0.17438]), np.array([-0.81301, 0.086, -0.77236]),
        np.array([-0.67901, 0.6566, -0.45438]), np.array([-0.02669, -0.13942, -0.98892]),
        np.array([1.01729, 0.83942, 1.15976]),
        np.array([1.92943, 1.10829, 0.87557]),
        np.array([0.58667, 0.97031, 1.31261]),
        np.array([1.57111, 0.74218, 0.68171]),
        1, 1, 1, 1, 2.0)
    e2 = get_erf_repulsion(
        np.array([0.13992, 0.37329, 0.33259]), np.array([0.64139, 1.73019, 0.13917]),
        np.array([0.44337, 1.28161, 0.3277]), np.array([1.24252, 1.27924, 1.45445]),
        np.array([0.02582, 0.94923, -0.17438]), np.array([-0.81301, 0.086, -0.77236]),
        np.array([-0.67901, 0.6566, -0.45438]), np.array([-0.02669, -0.13942, -0.98892]),
        np.array([1.01729, 0.83942, 1.15976]),
        np.array([1.92943, 1.10829, 0.87557]),
        np.array([0.58667, 0.97031, 1.31261]),
        np.array([1.57111, 0.74218, 0.68171]),
        1, 1, 1, 1, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([0.13992, 0.37329, 0.33259]), np.array([0.64139, 1.73019, 0.13917]),
        np.array([0.44337, 1.28161, 0.3277]), np.array([1.24252, 1.27924, 1.45445]),
        np.array([0.02582, 0.94923, -0.17438]), np.array([-0.81301, 0.086, -0.77236]),
        np.array([-0.67901, 0.6566, -0.45438]), np.array([-0.02669, -0.13942, -0.98892]),
        np.array([1.01729, 0.83942, 1.15976]),
        np.array([1.92943, 1.10829, 0.87557]),
        np.array([0.58667, 0.97031, 1.31261]),
        np.array([1.57111, 0.74218, 0.68171]),
        1, 1, 1, 1, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_0_2_1_0():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([1.36794, 1.14001, 1.97798]), np.array([1.68538, 0.75019, 0.72741]),
        np.array([1.55248, 0.78842, 1.84644]), np.array([1.73266, 0.46153, 0.63621]),
        np.array([0.05517, 0.27196, -0.98928]), np.array([-0.20526, 0.27314, -0.16208]),
        np.array([-0.00876, -0.47585, 0.88613]), np.array([0.75034, 0.54371, -0.1464]),
        np.array([0.50974]),
        np.array([1.36246, 0.58913, 0.73488, 0.53568, 1.11864, 1.80388]),
        np.array([1.62815, 0.58942, 1.52452]),
        np.array([0.66094]),
        0, 2, 1, 0, 2.0)
    e2 = get_erf_repulsion(
        np.array([1.36794, 1.14001, 1.97798]), np.array([1.68538, 0.75019, 0.72741]),
        np.array([1.55248, 0.78842, 1.84644]), np.array([1.73266, 0.46153, 0.63621]),
        np.array([0.05517, 0.27196, -0.98928]), np.array([-0.20526, 0.27314, -0.16208]),
        np.array([-0.00876, -0.47585, 0.88613]), np.array([0.75034, 0.54371, -0.1464]),
        np.array([0.50974]),
        np.array([1.36246, 0.58913, 0.73488, 0.53568, 1.11864, 1.80388]),
        np.array([1.62815, 0.58942, 1.52452]),
        np.array([0.66094]),
        0, 2, 1, 0, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([1.36794, 1.14001, 1.97798]), np.array([1.68538, 0.75019, 0.72741]),
        np.array([1.55248, 0.78842, 1.84644]), np.array([1.73266, 0.46153, 0.63621]),
        np.array([0.05517, 0.27196, -0.98928]), np.array([-0.20526, 0.27314, -0.16208]),
        np.array([-0.00876, -0.47585, 0.88613]), np.array([0.75034, 0.54371, -0.1464]),
        np.array([0.50974]),
        np.array([1.36246, 0.58913, 0.73488, 0.53568, 1.11864, 1.80388]),
        np.array([1.62815, 0.58942, 1.52452]),
        np.array([0.66094]),
        0, 2, 1, 0, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_0_2_2_3():
    epsilon = 0.000001
    e1 = get_erf_repulsion(
        np.array([0.96867, 0.41743, 1.03509]), np.array([1.84594, 0.83035, 1.20242]),
        np.array([0.94861, 0.47292, 0.38655]), np.array([1.3009, 1.10486, 1.4979]),
        np.array([0.10017, 0.21708, 0.08942]), np.array([-0.03049, 0.99486, -0.37959]),
        np.array([-0.7765, 0.53988, 0.25643]), np.array([0.60758, 0.85146, 0.15088]),
        np.array([1.14284]),
        np.array([1.39723, 1.77896, 0.72525, 0.99877, 1.5953, 0.69473]),
        np.array([0.56774, 1.69348, 1.8146, 0.85426, 1.35434, 1.87402]),
        np.array([0.99964, 1.45499, 1.35143, 1.9758, 0.58887, 1.40713, 0.55226,
                  1.44979, 0.57156, 0.71009]),
        0, 2, 2, 3, 2.0)
    e2 = get_erf_repulsion(
        np.array([0.96867, 0.41743, 1.03509]), np.array([1.84594, 0.83035, 1.20242]),
        np.array([0.94861, 0.47292, 0.38655]), np.array([1.3009, 1.10486, 1.4979]),
        np.array([0.10017, 0.21708, 0.08942]), np.array([-0.03049, 0.99486, -0.37959]),
        np.array([-0.7765, 0.53988, 0.25643]), np.array([0.60758, 0.85146, 0.15088]),
        np.array([1.14284]),
        np.array([1.39723, 1.77896, 0.72525, 0.99877, 1.5953, 0.69473]),
        np.array([0.56774, 1.69348, 1.8146, 0.85426, 1.35434, 1.87402]),
        np.array([0.99964, 1.45499, 1.35143, 1.9758, 0.58887, 1.40713, 0.55226,
                  1.44979, 0.57156, 0.71009]),
        0, 2, 2, 3, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([0.96867, 0.41743, 1.03509]), np.array([1.84594, 0.83035, 1.20242]),
        np.array([0.94861, 0.47292, 0.38655]), np.array([1.3009, 1.10486, 1.4979]),
        np.array([0.10017, 0.21708, 0.08942]), np.array([-0.03049, 0.99486, -0.37959]),
        np.array([-0.7765, 0.53988, 0.25643]), np.array([0.60758, 0.85146, 0.15088]),
        np.array([1.14284]),
        np.array([1.39723, 1.77896, 0.72525, 0.99877, 1.5953, 0.69473]),
        np.array([0.56774, 1.69348, 1.8146, 0.85426, 1.35434, 1.87402]),
        np.array([0.99964, 1.45499, 1.35143, 1.9758, 0.58887, 1.40713, 0.55226, 1.44979,
                  0.57156, 0.71009]),
        0, 2, 2, 3, egauss, 1.1283791670955126, 4.0)


def test_gauss_repulsion_4_3_2_1():
    epsilon = 0.0000001
    e1 = get_erf_repulsion(
        np.array([0.94212, 1.71823, 0.3309]), np.array([0.94854, 0.12816, 0.42016]),
        np.array([0.46046, 0.43321, 1.0587]), np.array([1.0089, 0.52286, 1.83539]),
        np.array([-0.48859, 0.6043, -0.57858]), np.array([0.74567, -0.82555, -0.30631]),
        np.array([-0.5679, -0.08725, 0.7623]), np.array([0.10338, 0.65407, -0.20172]),
        np.array([1.10904, 1.40637, 1.8707, 0.68295, 1.29692, 0.99892, 1.13936, 0.81258,
                  0.50325, 1.27698, 1.81192, 1.43415, 1.1686, 1.38063, 0.61592]),
        np.array([1.19368, 0.75291, 0.63535, 1.22654, 1.32848, 1.17482, 1.74897, 0.93964,
                  1.90303, 1.44528]),
        np.array([1.63343, 1.80498, 1.61313, 0.99992, 1.04505, 1.42297]),
        np.array([1.4825, 1.69421, 1.8635]),
        4, 3, 2, 1, 2.0)
    e2 = get_erf_repulsion(
        np.array([0.94212, 1.71823, 0.3309]), np.array([0.94854, 0.12816, 0.42016]),
        np.array([0.46046, 0.43321, 1.0587]), np.array([1.0089, 0.52286, 1.83539]),
        np.array([-0.48859, 0.6043, -0.57858]), np.array([0.74567, -0.82555, -0.30631]),
        np.array([-0.5679, -0.08725, 0.7623]), np.array([0.10338, 0.65407, -0.20172]),
        np.array([1.10904, 1.40637, 1.8707, 0.68295, 1.29692, 0.99892, 1.13936, 0.81258,
                  0.50325, 1.27698, 1.81192, 1.43415, 1.1686, 1.38063, 0.61592]),
        np.array([1.19368, 0.75291, 0.63535, 1.22654, 1.32848, 1.17482, 1.74897, 0.93964,
                  1.90303, 1.44528]),
        np.array([1.63343, 1.80498, 1.61313, 0.99992, 1.04505, 1.42297]),
        np.array([1.4825, 1.69421, 1.8635]),
        4, 3, 2, 1, 2.0 + epsilon)
    egauss = (e2 - e1) / epsilon
    check_gauss_repulsion(
        np.array([0.94212, 1.71823, 0.3309]), np.array([0.94854, 0.12816, 0.42016]),
        np.array([0.46046, 0.43321, 1.0587]), np.array([1.0089, 0.52286, 1.83539]),
        np.array([-0.48859, 0.6043, -0.57858]), np.array([0.74567, -0.82555, -0.30631]),
        np.array([-0.5679, -0.08725, 0.7623]), np.array([0.10338, 0.65407, -0.20172]),
        np.array([1.10904, 1.40637, 1.8707, 0.68295, 1.29692, 0.99892, 1.13936, 0.81258,
                  0.50325, 1.27698, 1.81192, 1.43415, 1.1686, 1.38063, 0.61592]),
        np.array([1.19368, 0.75291, 0.63535, 1.22654, 1.32848, 1.17482, 1.74897, 0.93964,
                  1.90303, 1.44528]),
        np.array([1.63343, 1.80498, 1.61313, 0.99992, 1.04505, 1.42297]),
        np.array([1.4825, 1.69421, 1.8635]),
        4, 3, 2, 1, egauss, 1.1283791670955126, 4.0)


def check_ralpha_repulsion(alphas0, alphas1, alphas2, alphas3, r0, r1, r2, r3, scales0,
                           scales1, scales2, scales3, shell_type0, shell_type1,
                           shell_type2, shell_type3, result0, alpha):
    """Compare output from HORTON Erf integrals with reference data.

    The reference data was generated with a Mathematica script of Julien Toulouse and
    Andreas Savin.

    Parameters
    ----------
    alpha0, alpha1, alpha2, alpha3 : float
        Exponents of the four primitive shells.
    r0, r1, r2, r3 : np.ndarray, shape=(3,), dtype=float
        Cartesian coordinates of the centers of the four primitive shells.
    scales0, scales1, scales2, scales3 : np.ndarray, dtype=float
        Normalization prefactors for the Gaussian shells.
    shell_type0, shell_type1, shell_type2, shell_type3 : int
        Shell types of the four primitive shells.
    result0 : np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis), dtype=float
        The expected result.
    alpha : float
        The interaction is r to the power alpha.
    """
    max_shell_type = 4
    max_nbasis = _get_shell_nbasis(max_shell_type)
    gb4i = _GB4RAlphaIntegralLibInt(max_shell_type, alpha)
    assert gb4i.max_nbasis == max_nbasis
    assert gb4i.nwork == max_nbasis ** 4
    assert gb4i.alpha == alpha

    nbasis0 = _get_shell_nbasis(shell_type0)
    nbasis1 = _get_shell_nbasis(shell_type1)
    nbasis2 = _get_shell_nbasis(shell_type2)
    nbasis3 = _get_shell_nbasis(shell_type3)
    assert result0.shape == (nbasis0, nbasis1, nbasis2, nbasis3)
    # Clear the working memory
    gb4i.reset(shell_type0, shell_type1, shell_type2, shell_type3, r0, r1, r2, r3)
    # Add a few contributions:
    for alpha0, alpha1, alpha2, alpha3 in zip(alphas0, alphas1, alphas2, alphas3):
        gb4i.add(1.0, alpha0, alpha1, alpha2, alpha3, scales0, scales1, scales2, scales3)
    result1 = gb4i.get_work(nbasis0, nbasis1, nbasis2, nbasis3)
    assert abs(result1 - result0).max() < 3e-7


def test_ralpha_simple0():
    check_ralpha_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[4.37335457]]]]), -1.)


def test_ralpha_simple1():
    check_ralpha_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[5.81367687]]]]), 2.)


def test_ralpha_simple2():
    check_ralpha_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[7.75156917]]]]), -2.)


def test_ralpha_repulsion_0_0_0_0_simple1():
    check_ralpha_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[2.20567322]]]]), -1.)


def test_ralpha_repulsion_0_0_0_0_simple2():
    check_ralpha_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0.57092, 0.29608, -0.758]), np.array([-0.70841, 0.22864, 0.79589]),
        np.array([0.83984, 0.65053, 0.36087]), np.array([-0.62267, -0.83676, -0.75233]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.19609589]]]]), -1.)


def test_ralpha_repulsion_0_0_0_0_simple3():
    check_ralpha_repulsion(
        np.array([0.57283]), np.array([1.74713]),
        np.array([0.21032]), np.array([1.60538]),
        np.array([0.82197, 0.73226, -0.98154]), np.array([0.57466, 0.17815, -0.25519]),
        np.array([0.00425, -0.33757, 0.08556]), np.array([-0.38717, 0.66721, 0.40838]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.92553047]]]]), -1.)


def test_ralpha_repulsion_0_0_0_0_simple4():
    check_ralpha_repulsion(
        np.array([1.35491]), np.array([0.9714]),
        np.array([1.95585]), np.array([1.77853]),
        np.array([0.37263, -0.87382, 0.28078]), np.array([-0.08946, -0.52616, 0.69184]),
        np.array([-0.35128, 0.07017, 0.08193]), np.array([0.14543, -0.29499, -0.09769]),
        np.array([1.61086]),
        np.array([1.19397]),
        np.array([1.8119]),
        np.array([1.55646]),
        0, 0, 0, 0,
        np.array([[[[1.65373353]]]]), -1.)


def test_ralpha_repulsion_0_0_0_0():
    check_ralpha_repulsion(
        np.array([1.63216, 1.25493, 1.46134, 0.48024]),
        np.array([1.72365, 1.59905, 0.10447, 1.28324]),
        np.array([1.4105, 0.27134, 1.51238, 0.7518]),
        np.array([1.38488, 0.97611, 0.34149, 0.4326]),
        np.array([0.61356, -0.85284, -0.37151]), np.array([-0.63238, -0.81396, 0.40314]),
        np.array([0.29559, 0.60342, 0.18878]), np.array([-0.6893, 0.09175, -0.97283]),
        np.array([0.84965]),
        np.array([1.49169]),
        np.array([1.11046]),
        np.array([0.6665]),
        0, 0, 0, 0,
        np.array([[[[2.97326773]]]]), -1.)


def test_ralpha_repulsion_0_0_0_1():
    check_ralpha_repulsion(
        np.array([0.74579, 0.93686, 0.39742]), np.array([1.01349, 1.46072, 0.22295]),
        np.array([1.90756, 0.52423, 1.35586]), np.array([0.9655, 0.73539, 0.51017]),
        np.array([0.55177, 0.11232, -0.95152]), np.array([0.79941, 0.80782, 0.02287]),
        np.array([-0.52471, 0.59124, 0.434]), np.array([0.40758, 0.96818, 0.59852]),
        np.array([1.38989]),
        np.array([1.20619]),
        np.array([1.25917]),
        np.array([0.70246, 1.69253, 1.5632]),
        0, 0, 0, 1,
        np.array([[[[0.03331693, -3.27978915, -5.8871596]]]]), -1.)


def test_ralpha_repulsion_0_0_1_0():
    check_ralpha_repulsion(
        np.array([0.88609, 0.76883, 0.56082]), np.array([1.29216, 0.28671, 1.25389]),
        np.array([1.36987, 0.90792, 0.30511]), np.array([0.57079, 1.98163, 0.66835]),
        np.array([0.7706, 0.99091, -0.21592]), np.array([-0.00566, -0.37522, -0.3936]),
        np.array([0.1527, -0.95347, 0.16682]), np.array([-0.75347, -0.6388, -0.81567]),
        np.array([0.72776]),
        np.array([1.08088]),
        np.array([0.53874, 1.37722, 1.16945]),
        np.array([1.4472]),
        0, 0, 1, 0,
        np.array([[[[0.5110845], [5.86518324], [-1.47878266]]]]), -1.0)


def test_ralpha_repulsion_0_0_1_1():
    check_ralpha_repulsion(
        np.array([0.94138, 0.23708, 1.33464]), np.array([1.89753, 0.54214, 0.80346]),
        np.array([1.04131, 1.6925, 0.81454]), np.array([1.06467, 0.55116, 1.21121]),
        np.array([0.6941, 0.3354, -0.49162]), np.array([0.68756, 0.49975, -0.69756]),
        np.array([0.60432, -0.01449, -0.26057]), np.array([0.35763, -0.04674, -0.78137]),
        np.array([0.75847]),
        np.array([0.57683]),
        np.array([1.61747, 0.59289, 0.93361]),
        np.array([1.38523, 1.77715, 0.8249]),
        0, 0, 1, 1,
        np.array([[[[0.57999607, 0.04732015, 0.00079488],
                    [0.09549513, 0.42707461, 0.03630467],
                    [-0.15902635, -0.25704193, 0.12295133]]]]),
        -1.0)


def test_ralpha_repulsion_0_1_0_0():
    check_ralpha_repulsion(
        np.array([0.11308, 0.49861, 1.12215]), np.array([0.6186, 1.93501, 1.72751]),
        np.array([0.4644, 0.61371, 1.99408]), np.array([1.98686, 0.49338, 0.88466]),
        np.array([0.31794, 0.18412, 0.89808]), np.array([0.35463, 0.17042, 0.0682]),
        np.array([0.51676, -0.86674, -0.32785]), np.array([-0.03453, -0.05741, -0.86135]),
        np.array([1.84487]),
        np.array([1.17293, 1.02836, 0.50605]),
        np.array([0.54734]),
        np.array([1.55774]),
        0, 1, 0, 0,
        np.array([[[[-2.98984233]], [[-2.16665085]], [[-3.19087757]]]]), -1.0)


def test_ralpha_repulsion_0_1_0_1():
    check_ralpha_repulsion(
        np.array([0.95345, 1.7616, 0.62144]), np.array([0.60537, 0.78954, 0.17662]),
        np.array([1.39946, 1.03161, 1.42837]), np.array([1.05228, 1.80211, 1.37614]),
        np.array([0.18086, -0.0927, -0.36495]), np.array([0.48062, -0.97782, -0.05878]),
        np.array([-0.55927, -0.95238, 0.33122]), np.array([0.17856, 0.06077, 0.62697]),
        np.array([0.9876]),
        np.array([1.39633, 1.30787, 1.80682]),
        np.array([0.93201]),
        np.array([1.21516, 1.84023, 1.59345]),
        0, 1, 0, 1,
        np.array([[[[1.11620596, 0.60061237, 0.36843148]],
                   [[-0.05340867, 0.33119515, -0.70418275]],
                   [[-0.04504112, -1.01394262, 1.17313632]]]]),
        -1.0)


def test_ralpha_repulsion_0_1_1_1():
    check_ralpha_repulsion(
        np.array([1.60961, 1.48434, 1.09022]), np.array([1.49016, 0.78972, 1.01383]),
        np.array([1.357, 1.6929, 1.46297]), np.array([1.3126, 1.39773, 0.3295]),
        np.array([-0.74441, 0.13168, 0.17287]), np.array([-0.73242, 0.73598, -0.07688]),
        np.array([0.06303, 0.61361, 0.92689]), np.array([0.31395, 0.00081, -0.13425]),
        np.array([1.92653]),
        np.array([0.84324, 1.68215, 0.64055]),
        np.array([1.62317, 1.94784, 1.54325]),
        np.array([0.67873, 0.76053, 0.57816]),
        0, 1, 1, 1,
        np.array([[[[-0.06633908, -0.13761956, -0.03005655],
                    [-0.023407, -0.07813472, -0.03489736],
                    [-0.02263273, -0.20143856, -0.03550443]],
                   [[-0.40044718, -0.35436776, 0.07827812],
                    [-0.39382673, -0.18295174, 0.10845718],
                    [-0.37310311, -0.34400264, 0.05152883]],
                   [[0.07743294, -0.04648822, -0.2043075],
                    [0.03540926, -0.00400861, -0.13446393],
                    [0.02364929, -0.01807209, -0.18079094]]]]),
        -1.0)


def test_ralpha_repulsion_1_0_0_1():
    check_ralpha_repulsion(
        np.array([0.39834, 1.4798, 1.80662]), np.array([1.9623, 0.88607, 0.93517]),
        np.array([0.46864, 1.1317, 0.67625]), np.array([1.52214, 0.93879, 0.71425]),
        np.array([-0.04796, 0.70504, 0.36481]), np.array([0.40599, 0.97607, 0.64758]),
        np.array([0.66271, -0.64123, -0.17474]), np.array([-0.60087, 0.25093, 0.32664]),
        np.array([0.68301, 1.18047, 1.44482]),
        np.array([0.97181]),
        np.array([1.18315]),
        np.array([0.79184, 1.41932, 1.32812]),
        1, 0, 0, 1,
        np.array([[[[0.16173756, 0.14265052, 0.05405344]]],
                  [[[-0.431925, -0.37295006, -0.1782411]]],
                  [[[-0.17915755, -0.20235955, 0.03526912]]]]),
        -1.0)


def test_ralpha_repulsion_1_1_1_1():
    check_ralpha_repulsion(
        np.array([0.13992, 0.37329, 0.33259]), np.array([0.64139, 1.73019, 0.13917]),
        np.array([0.44337, 1.28161, 0.3277]), np.array([1.24252, 1.27924, 1.45445]),
        np.array([0.02582, 0.94923, -0.17438]), np.array([-0.81301, 0.086, -0.77236]),
        np.array([-0.67901, 0.6566, -0.45438]), np.array([-0.02669, -0.13942, -0.98892]),
        np.array([1.01729, 0.83942, 1.15976]),
        np.array([1.92943, 1.10829, 0.87557]),
        np.array([0.58667, 0.97031, 1.31261]),
        np.array([1.57111, 0.74218, 0.68171]),
        1, 1, 1, 1,
        np.array([[[[5.38092832, 0.67101024, 0.50643354],
                    [-0.36637823, -0.17128347, 0.00749151],
                    [-0.47015285, -0.00846274, -0.23514519]],
                   [[0.31412053, 1.85552661, -0.05096966],
                    [0.5668773, 0.04019152, -0.05803149],
                    [-0.02195855, 0.00256108, 0.03373068]],
                   [[0.26139911, -0.05908764, 1.34729127],
                    [0.03563575, 0.02599451, 0.0669569],
                    [0.6249628, -0.09012696, -0.02559206]]],
                  [[[-1.30079959, 0.06525516, -0.24130176],
                    [7.90805546, 0.5029288, 1.03164863],
                    [0.22531828, -0.01518479, -0.63472654]],
                   [[0.07758755, -0.30344079, 0.03679751],
                    [0.88274549, 3.43263474, -0.20761467],
                    [0.09249023, 0.10854722, 0.15741632]],
                   [[0.0082139, -0.00382022, -0.24202072],
                    [0.44155444, -0.06437548, 2.40552259],
                    [0.29276089, 0.01725224, 0.05956368]]],
                  [[[-1.45339037, -0.37266055, 0.25844897],
                    [0.41152374, -0.40525461, -0.16607501],
                    [14.23224926, 2.34068558, 0.65653732]],
                   [[-0.00776144, -0.38261119, -0.0073076],
                    [0.28311943, 0.14089539, 0.08426703],
                    [0.91304633, 5.92042353, -0.12886949]],
                   [[0.09807363, 0.06281554, -0.25920407],
                    [0.15636252, 0.10752926, 0.14182457],
                    [1.2142302, -0.38098265, 4.57694241]]]]),
        -1.0)


def test_ralpha_repulsion_0_2_1_0():
    check_ralpha_repulsion(
        np.array([1.36794, 1.14001, 1.97798]), np.array([1.68538, 0.75019, 0.72741]),
        np.array([1.55248, 0.78842, 1.84644]), np.array([1.73266, 0.46153, 0.63621]),
        np.array([0.05517, 0.27196, -0.98928]), np.array([-0.20526, 0.27314, -0.16208]),
        np.array([-0.00876, -0.47585, 0.88613]), np.array([0.75034, 0.54371, -0.1464]),
        np.array([0.50974]),
        np.array([1.36246, 0.58913, 0.73488, 0.53568, 1.11864, 1.80388]),
        np.array([1.62815, 0.58942, 1.52452]),
        np.array([0.66094]),
        0, 2, 1, 0,
        np.array([[[[0.03940319], [0.05597157], [-0.32990373]],
                   [[-0.00066587], [0.00221213], [-0.00319745]],
                   [[-0.00035194], [-0.00011777], [0.0063613]],
                   [[0.00478058], [0.01592957], [-0.09687372]],
                   [[0.00002574], [-0.00009517], [-0.00166564]],
                   [[0.01578456], [0.05420504], [-0.32175899]]]]),
        -1.0)


def test_ralpha_repulsion_0_2_2_3():
    fn = 'electron_repulsion_0_2_2_3_json'
    result0 = load_json(fn)
    check_ralpha_repulsion(
        np.array([0.96867, 0.41743, 1.03509]), np.array([1.84594, 0.83035, 1.20242]),
        np.array([0.94861, 0.47292, 0.38655]), np.array([1.3009, 1.10486, 1.4979]),
        np.array([0.10017, 0.21708, 0.08942]), np.array([-0.03049, 0.99486, -0.37959]),
        np.array([-0.7765, 0.53988, 0.25643]), np.array([0.60758, 0.85146, 0.15088]),
        np.array([1.14284]),
        np.array([1.39723, 1.77896, 0.72525, 0.99877, 1.5953, 0.69473]),
        np.array([0.56774, 1.69348, 1.8146, 0.85426, 1.35434, 1.87402]),
        np.array([0.99964, 1.45499, 1.35143, 1.9758, 0.58887, 1.40713, 0.55226, 1.44979,
                  0.57156, 0.71009]),
        0, 2, 2, 3,
        result0, -1.0)


def test_ralpha_repulsion_4_3_2_1():
    fn = 'electron_repulsion_4_3_2_1_json'
    result0 = load_json(fn)
    check_ralpha_repulsion(
        np.array([0.94212, 1.71823, 0.3309]), np.array([0.94854, 0.12816, 0.42016]),
        np.array([0.46046, 0.43321, 1.0587]), np.array([1.0089, 0.52286, 1.83539]),
        np.array([-0.48859, 0.6043, -0.57858]), np.array([0.74567, -0.82555, -0.30631]),
        np.array([-0.5679, -0.08725, 0.7623]), np.array([0.10338, 0.65407, -0.20172]),
        np.array([1.10904, 1.40637, 1.8707, 0.68295, 1.29692, 0.99892, 1.13936, 0.81258,
                  0.50325, 1.27698, 1.81192, 1.43415, 1.1686, 1.38063, 0.61592]),
        np.array([1.19368, 0.75291, 0.63535, 1.22654, 1.32848, 1.17482, 1.74897, 0.93964,
                  1.90303, 1.44528]),
        np.array([1.63343, 1.80498, 1.61313, 0.99992, 1.04505, 1.42297]),
        np.array([1.4825, 1.69421, 1.8635]),
        4, 3, 2, 1,
        result0, -1.0)


def get_gauss_repulsion(alphas0, alphas1, alphas2, alphas3, r0, r1, r2, r3,
                        scales0, scales1, scales2, scales3,
                        shell_type0, shell_type1, shell_type2, shell_type3, c, alpha):
    """Get the short-range damped Erf integrals for a primitive shell.

    Parameters
    ----------
    alpha0, alpha1, alpha2, alpha3 : float
        Exponents of the four primitive shells.
    r0, r1, r2, r3 : np.ndarray, shape=(3,), dtype=float
        Cartesian coordinates of the centers of the four primitive shells.
    scales0, scales1, scales2, scales3 : float
        Normalization prefactors for the Gaussian shells.
    shell_type0, shell_type1, shell_type2, shell_type3 : int
        Shell types of the four primitive shells.
    c : float
        The coefficient of the Gaussian function.
    alpha : float
        The exponent of the Gaussian function.

    """
    max_shell_type = 4
    gb4i = _GB4GaussIntegralLibInt(max_shell_type, c, alpha)

    nbasis0 = _get_shell_nbasis(shell_type0)
    nbasis1 = _get_shell_nbasis(shell_type1)
    nbasis2 = _get_shell_nbasis(shell_type2)
    nbasis3 = _get_shell_nbasis(shell_type3)
    # Clear the working memory
    gb4i.reset(shell_type0, shell_type1, shell_type2, shell_type3, r0, r1, r2, r3)
    # Add a few cobtributions:
    for alpha0, alpha1, alpha2, alpha3 in zip(alphas0, alphas1, alphas2, alphas3):
        gb4i.add(1.0, alpha0, alpha1, alpha2, alpha3, scales0, scales1, scales2, scales3)
    return gb4i.get_work(nbasis0, nbasis1, nbasis2, nbasis3)


def check_delta_repulsion(alphas0, alphas1, alphas2, alphas3, r0, r1, r2, r3, scales0,
                          scales1, scales2, scales3, shell_type0, shell_type1,
                          shell_type2, shell_type3, result0):
    """Compare output from HORTON Delta 4-center integrals with reference data.

    The reference data was generated using very pointy Gauss 4-center integrals.

    Parameters
    ----------
    alpha0, alpha1, alpha2, alpha3 : float
        Exponents of the four primitive shells.
    r0, r1, r2, r3 : np.ndarray, shape=(3,), dtype=float
        Cartesian coordinates of the centers of the four primitive shells.
    scales0, scales1, scales2, scales3 : float
        Normalization prefactors for the Gaussian shells.
    shell_type0, shell_type1, shell_type2, shell_type3 : int
        Shell types of the four primitive shells.
    result0 : np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis), dtype=float
        The expected result.

    """
    max_shell_type = 4
    max_nbasis = _get_shell_nbasis(max_shell_type)
    gb4i = _GB4DeltaIntegralLibInt(max_shell_type)
    assert gb4i.max_nbasis == max_nbasis
    assert gb4i.nwork == max_nbasis ** 4

    nbasis0 = _get_shell_nbasis(shell_type0)
    nbasis1 = _get_shell_nbasis(shell_type1)
    nbasis2 = _get_shell_nbasis(shell_type2)
    nbasis3 = _get_shell_nbasis(shell_type3)
    assert result0.shape == (nbasis0, nbasis1, nbasis2, nbasis3)
    # Clear the working memory
    gb4i.reset(shell_type0, shell_type1, shell_type2, shell_type3, r0, r1, r2, r3)
    # Add a few cobtributions:
    for alpha0, alpha1, alpha2, alpha3 in zip(alphas0, alphas1, alphas2, alphas3):
        gb4i.add(1.0, alpha0, alpha1, alpha2, alpha3, scales0, scales1, scales2, scales3)
    result1 = gb4i.get_work(nbasis0, nbasis1, nbasis2, nbasis3)
    print("number ", result1)
    assert abs(result1 - result0).max() < 3e-7


def test_delta_simple0():
    check_delta_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.6960409996]]]]))


def test_delta_repulsion_0_0_0_0_simple1():
    check_delta_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.03465384]]]]))


def test_delta_repulsion_0_0_0_0_simple2():
    check_delta_repulsion(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0.57092, 0.29608, -0.758]), np.array([-0.70841, 0.22864, 0.79589]),
        np.array([0.83984, 0.65053, 0.36087]), np.array([-0.62267, -0.83676, -0.75233]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.00456547]]]]))


def test_delta_repulsion_0_0_0_0_simple3():
    check_delta_repulsion(
        np.array([0.57283]), np.array([1.74713]),
        np.array([0.21032]), np.array([1.60538]),
        np.array([0.82197, 0.73226, -0.98154]), np.array([0.57466, 0.17815, -0.25519]),
        np.array([0.00425, -0.33757, 0.08556]), np.array([-0.38717, 0.66721, 0.40838]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[[[0.06552595]]]]))


def test_delta_repulsion_0_0_0_0_simple4():
    check_delta_repulsion(
        np.array([1.35491]), np.array([0.9714]),
        np.array([1.95585]), np.array([1.77853]),
        np.array([0.37263, -0.87382, 0.28078]), np.array([-0.08946, -0.52616, 0.69184]),
        np.array([-0.35128, 0.07017, 0.08193]), np.array([0.14543, -0.29499, -0.09769]),
        np.array([1.61086]),
        np.array([1.19397]),
        np.array([1.8119]),
        np.array([1.55646]),
        0, 0, 0, 0,
        np.array([[[[0.3883875489]]]]))


# Skip test. Probably wrong.
# def test_delta_repulsion_0_0_0_1():
#     sigma = 0.0001
#     c = 1.0 / (np.sqrt(2.0 * np.pi) * sigma)
#     alpha = 1.0 / (2.0 * sigma * sigma)
#     e1 = get_gauss_repulsion(
#         np.array([0.74579, 0.93686, 0.39742]), np.array([1.01349, 1.46072, 0.22295]),
#         np.array([1.90756, 0.52423, 1.35586]), np.array([0.9655, 0.73539, 0.51017]),
#         np.array([0.55177, 0.11232, -0.95152]), np.array([0.79941, 0.80782, 0.02287]),
#         np.array([-0.52471, 0.59124, 0.434]), np.array([0.40758, 0.96818, 0.59852]),
#         np.array([1.38989]),
#         np.array([1.20619]),
#         np.array([1.25917]),
#         np.array([0.70246, 1.69253, 1.5632]),
#         0, 0, 0, 1, c, alpha)
#     check_delta_repulsion(
#         np.array([0.74579, 0.93686, 0.39742]), np.array([1.01349, 1.46072, 0.22295]),
#         np.array([1.90756, 0.52423, 1.35586]), np.array([0.9655, 0.73539, 0.51017]),
#         np.array([0.55177, 0.11232, -0.95152]), np.array([0.79941, 0.80782, 0.02287]),
#         np.array([-0.52471, 0.59124, 0.434]), np.array([0.40758, 0.96818, 0.59852]),
#         np.array([1.38989]),
#         np.array([1.20619]),
#         np.array([1.25917]),
#         np.array([0.70246, 1.69253, 1.5632]),
#         0, 0, 0, 1, e1)


def check_intracule_integrals(alphas0, alphas1, alphas2, alphas3, r0, r1, r2, r3, scales0,
                              scales1, scales2, scales3, shell_type0, shell_type1,
                              shell_type2, shell_type3, point, result0):
    """Compare output from HORTON Delta 4-center integrals with reference data.

    The reference data was generated using very pointy Gauss 4-center integrals.

    Parameters
    ----------
    alpha0, alpha1, alpha2, alpha3 : float
        Exponents of the four primitive shells.
    r0, r1, r2, r3 : np.ndarray, shape=(3,), dtype=float
        Cartesian coordinates of the centers of the four primitive shells.
    scales0, scales1, scales2, scales3 : float
        Normalization prefactors for the Gaussian shells.
    shell_type0, shell_type1, shell_type2, shell_type3 : int
        Shell types of the four primitive shells.
    result0 : np.ndarray, shape=(nbasis, nbasis, nbasis, nbasis), dtype=float
        The expected result.

    """
    max_shell_type = 4
    max_nbasis = _get_shell_nbasis(max_shell_type)
    gb4i = _GB4IntraDensIntegralLibInt(max_shell_type, point)
    assert gb4i.max_nbasis == max_nbasis
    assert gb4i.nwork == max_nbasis ** 4

    nbasis0 = _get_shell_nbasis(shell_type0)
    nbasis1 = _get_shell_nbasis(shell_type1)
    nbasis2 = _get_shell_nbasis(shell_type2)
    nbasis3 = _get_shell_nbasis(shell_type3)
    assert result0.shape == (nbasis0, nbasis1, nbasis2, nbasis3)
    # Clear the working memory
    gb4i.reset(shell_type0, shell_type1, shell_type2, shell_type3, r0, r1, r2, r3)
    # Add a few cobtributions:
    for alpha0, alpha1, alpha2, alpha3 in zip(alphas0, alphas1, alphas2, alphas3):
        gb4i.add(1.0, alpha0, alpha1, alpha2, alpha3, scales0, scales1, scales2, scales3)
    result1 = gb4i.get_work(nbasis0, nbasis1, nbasis2, nbasis3)
    print("inside result ", result1)
    assert abs(result1 - result0).max() < 3e-7


def test_intracule_simple0():
    check_intracule_integrals(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([0., 0., 0.]), np.array([0., 0., 0.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[0., 1., 0.]]),
        np.array([[[[0.25605917396]]]]))


def test_intracule_repulsion_0_0_0_0_simple1():
    check_intracule_integrals(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([0., 0., 0.]), np.array([1., 1., 1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[0., 1., 0.]]),
        np.array([[[[0.00172531314]]]]))


def test_intracule_repulsion_0_0_0_0_simple2():
    check_intracule_integrals(
        np.array([1.]), np.array([1.]),
        np.array([1.]), np.array([1.]),
        np.array([0.57092, 0.29608, -0.758]), np.array([-0.70841, 0.22864, 0.79589]),
        np.array([0.83984, 0.65053, 0.36087]), np.array([-0.62267, -0.83676, -0.75233]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[0., 1., 0.]]),
        np.array([[[[0.0079506238]]]]))


def test_intracule_repulsion_0_0_0_0_simple3():
    check_intracule_integrals(
        np.array([0.57283]), np.array([1.74713]),
        np.array([0.21032]), np.array([1.60538]),
        np.array([0.82197, 0.73226, -0.98154]), np.array([0.57466, 0.17815, -0.25519]),
        np.array([0.00425, -0.33757, 0.08556]), np.array([-0.38717, 0.66721, 0.40838]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        np.array([1.]),
        0, 0, 0, 0,
        np.array([[0., 1., 0.]]),
        np.array([[[[0.036197920998]]]]))


def test_intracule_repulsion_0_0_0_0_simple4():
    check_intracule_integrals(
        np.array([1.35491]), np.array([0.9714]),
        np.array([1.95585]), np.array([1.77853]),
        np.array([0.37263, -0.87382, 0.28078]), np.array([-0.08946, -0.52616, 0.69184]),
        np.array([-0.35128, 0.07017, 0.08193]), np.array([0.14543, -0.29499, -0.09769]),
        np.array([1.61086]),
        np.array([1.19397]),
        np.array([1.8119]),
        np.array([1.55646]),
        0, 0, 0, 0,
        np.array([[0., 1., 0.]]),
        np.array([[[[0.10370627969]]]]))


def test_delta_repulsion_int0():
    fn = 'water_ccpvdz_pure_hf_g03_fchk'
    obasis = load_obasis(fn)
    obasis.compute_delta_repulsion(shift=np.zeros([4, 3]))


def test_delta_repulsion_int1():
    coords = np.array([[0., 0., 1.], [0., 1., 0.], [0., 1., 1.], [0., 0., 0.]])
    numbers = np.ones([4], dtype=int)
    basis = "sto-3g"
    gobasis_ref = get_gobasis(coords, numbers, basis)

    coords1 = np.array([[0, 0, 2], [0, 2, 0], [0, 2, 2], [1, 1, 1]], dtype=float)
    shift = np.array([[0, 0, -1], [0, -1, 0], [0, -1, -1], [-1, -1, -1, ]], dtype=float)
    gobasis = get_gobasis(coords1, numbers, basis)

    assert np.allclose(gobasis.compute_delta_repulsion(shift=shift),
                       gobasis_ref.compute_delta_repulsion())
