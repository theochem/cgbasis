# -*- coding: utf-8 -*-
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


from gbasis.test.cext import _fac, _fac2, _jfac, _cit, _binom


def test_fac():
    assert _fac(-20) == 1
    assert _fac(0) == 1
    assert _fac(1) == 1
    assert _fac(2) == 2
    assert _fac(3) == 6
    assert _fac(4) == 24
    assert _fac(5) == 120
    assert _fac(10) == 3628800


def test_fac2():
    assert _fac2(-20) == 1
    assert _fac2(0) == 1
    assert _fac2(1) == 1
    assert _fac2(2) == 2
    assert _fac2(3) == 3
    assert _fac2(4) == 8
    assert _fac2(5) == 15


def test_binom():
    assert _binom(1, 1) == 1
    assert _binom(5, 3) == 10
    assert _binom(3, 2) == 3
    assert _binom(10, 4) == 210
    assert _binom(18, 14) == 3060
    assert _binom(5, 1) == 5
    assert _binom(5, 0) == 1
    assert _binom(0, 0) == 1
    assert _binom(5, 5) == 1


def test_cit():
    assert abs(_cit(1, 2.49128, 1) - 1.6608533333) < 1e-7
    assert abs(_cit(2, 2.49128, 2) - 0.8275301384) < 1e-7
    assert abs(_cit(4, 0.0213, 2) - 3.20063e-7) < 1e-9
    assert abs(_cit(2, 1.5134, 4) - 0.69944513718) < 1e-7


def test_jfac():
    assert _jfac(2, 3) == 0
    assert _jfac(3, 1) == 3
    assert _jfac(3, 3) == 6
    assert _jfac(5, 4) == 120
    assert _jfac(10, 3) == 720
    assert _jfac(10, 4) == 5040
