// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2017 The HORTON Development Team
//
// This file is part of HORTON.
//
// HORTON is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// HORTON is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--

/**
 * @file cartpure.h
 * @brief Conversion of Cartesian to Pure Gaussian functions
 */

#ifndef GBASIS_CARTPURE_H_
#define GBASIS_CARTPURE_H_

/** @brief
    Convert results for a Cartesian basis into results for a Pure basis, just
    by taking pre-programmed linear combinations.

    @param work_cart
        The input array with the Cartesian results in contiguous storage

    @param work_pure
        The output array for the Pure results, also in contiguous storage

    @param shell_type
        The Cartesian shell type that must be converted, i.e. >= 0

    @param nant, npost
        These two integer parameters describe how each record of Cartesian
        and Pure results is stored. A record is defined as a set of results
        for a single shell. The convention is as follows:

            work_cart[(ca*ncart + icart)*npost + cp]
            work_pure[(ca*npure + ipure)*npost + cp]

        where:

        * ca = anterior counter for the records (0..nant-1)
        * cp = anterior counter for the records (0..npost-1)
        * ncart = the number of Cartesian basis functions in this shell
        * npure = the number of Pure basis functions in this shell
        * icart = a counter for the Cartesian basis function in this shell
        * ipure = a counter for the Pure basis function in this shell
        * nant*npost = the number of records
        * npost = the stride between two values in the same record
*/
void cart_to_pure_low(double *work_cart, double *work_pure, long shell_type,
                      long nant, long npost);

#endif  // GBASIS_CARTPURE_H_
