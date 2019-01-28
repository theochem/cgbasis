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


from gbasis.pxds cimport c_gbasis

cdef extern from "grids/iter_gb1.h":
    cdef cppclass IterGB1:
        IterGB1(c_gbasis.GBasis* gbasis)

        bint inc_shell()
        void update_shell()
        bint inc_prim()
        void update_prim()
        void store(double* work, double* output, long dim)

        # 'public' iterator fields
        long shell_type0
        double con_coeff, alpha0
        double* r0
        double* scales0
        long ibasis0

        # 'private' iterator fields
        long ishell0
        long nprim0, iprim0, oprim0
