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
 * @file iter_pow1.h
 * @brief Iterators over Cartesian polynomials in one shell
 */

#ifndef GBASIS_ITER_POW1_H
#define GBASIS_ITER_POW1_H

///Increment to next basis function, within the same angular momentum.
int iter_pow1_inc(long *n);

/**
 * Iterator over 1 cartesian polynomial
 */
class IterPow1 {
 private:
  /// shell type
  long shell_type0;
 public:
  /// reinitialize on a new shell
  void reset(long shell_type0);
  /// increment to next basis function, first within the shell, then to the next shell.
  int inc();
  /// 3D cartesian shell type
  long n0[3];
  /// basis counter
  long ibasis0;
};

#endif
