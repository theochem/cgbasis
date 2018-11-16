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


from .. cimport c_gbasis

cdef extern from "ones/iter_gb2.h":
    cdef cppclass IterGB2:
        IterGB2(c_gbasis.GBasis* gbasis)

        bint inc_shell()
        void update_shell()
        bint inc_prim()
        void update_prim()
        void store(double* work, double* output)

        # 'public' iterator fields
        long shell_type0, shell_type1
        double con_coeff, alpha0, alpha1
        double* r0
        double* r1
        double* scales0
        double* scales1
        long ibasis0, ibasis1

        # 'private' iterator fields
        long ishell0, ishell1
        long nprim0, nprim1, iprim0, iprim1, oprim0, oprim1
