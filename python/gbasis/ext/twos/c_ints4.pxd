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


cdef extern from "twos/ints4.h":
    cdef cppclass GB4Integral:
        long get_nwork()
        long get_max_shell_type()
        long get_max_nbasis()
        void reset(long shell_type0, long shell_type1, long shell_type2, long shell_type3, double* r0, double* r1, double* r2, double* r3) except +
        void add(double coeff, double alpha0, double alpha1, double alpha2, double alpha3, double* scales0, double* scales1, double* scales2, double* scales3)
        void cart_to_pure() except +

        long get_shell_type0()
        long get_shell_type1()
        long get_shell_type2()
        long get_shell_type3()
        double* get_work()

    cdef cppclass GB4ElectronRepulsionIntegralLibInt:
        GB4ElectronRepulsionIntegralLibInt(long max_shell_type) except +

    cdef cppclass GB4ErfIntegralLibInt:
        GB4ErfIntegralLibInt(long max_shell_type, double mu) except +
        double get_mu()

    cdef cppclass GB4GaussIntegralLibInt:
        GB4GaussIntegralLibInt(long max_shell_type, double c, double alpha) except +
        double get_c()
        double get_alpha()

    cdef cppclass GB4RAlphaIntegralLibInt:
        GB4RAlphaIntegralLibInt(long max_shell_type, double alpha) except +
        double get_alpha()

    cdef cppclass GB4DeltaIntegralLibInt:
        GB4DeltaIntegralLibInt(long max_shell_type) except +

    cdef cppclass GB4IntraDensIntegralLibInt:
        GB4IntraDensIntegralLibInt(long max_shell_type, double* point) except +


cdef extern from "libint2.h":
    void libint2_static_init()
    void libint2_static_cleanup()
