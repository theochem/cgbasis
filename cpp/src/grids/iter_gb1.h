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
 * @file iter_gb.h
 * @brief Iterators over Gaussian basis functions
 */

#ifndef GBASIS_ITER_GB1_H_
#define GBASIS_ITER_GB1_H_

#include "../gbasis.h"
//! A 1 basis function iterator. Not used very much except to define basic structure.
class IterGB1 {
 private:
  // input fields
  const GBasis *gbasis;
  const long *basis_offsets;

 public:
  explicit IterGB1(GBasis *gbasis);
  IterGB1(const IterGB1 &other) = delete;

  //! Increments the shell index
  int inc_shell();

  void update_shell();

  int inc_prim();

  void update_prim();

  void store(const double *work, double *output, long dim);

  // 'public' iterator fields
  long shell_type0;
  double con_coeff, alpha0;
  const double *r0;  /// The current center
  const double *scales0;  /// Normalization constants
  long ibasis0;  /// Basis function counters (for output storage)

  // 'private' iterator fields
  long ishell0;  /// Shell counters
  long nprim0, iprim0, oprim0;  /// Primitive counters
};

#endif  // GBASIS_ITER_GB1_H_
