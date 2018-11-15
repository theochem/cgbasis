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

#ifndef GBASIS_ITER_GB4_H_
#define GBASIS_ITER_GB4_H_

#include "../gbasis.h"

//! Iterates over 4 basis functions in order. Used for 2-electron integrals.
class IterGB4 {
 private:
  // input fields
  const GBasis *gbasis;
  const long *basis_offsets;

 public:
  explicit IterGB4(GBasis *gbasis);
  IterGB4(const IterGB4 &other) = delete;

  /// Increment the shell counters (oprims, ishells) Takes into account 8-fold symmetry.
  int inc_shell();

  /// Update fields that depend on primitive counters (centres, shells types, nprim) from the gbasis instance.
  void update_shell();

  /// Increment the primitive counters (iprims)
  int inc_prim();

  /// Update fields that depend on primitive counters (alphas, contraction coeffs, and scaling factors) from gbasis.
  void update_prim();
  /**
   * Write the current shell, stored in work, into the output array. Only works for dense storage at the moment.
   *
   * @param work array of pure basis functions for the current shell.
   * @param output array storing the entire integral.
   */
  void store(const double *work, double *output);

  // 'public' iterator fields

  /// shell types
  long shell_type0, shell_type1, shell_type2, shell_type3;
  /// Contraction coefficient
  double con_coeff;
  /// Exponents
  double alpha0, alpha1, alpha2, alpha3;
  /// current gaussian centers
  const double *r0;
  const double *r1;
  const double *r2;
  const double *r3;
  /// normalization constants
  const double *scales0;
  const double *scales1;
  const double *scales2;
  const double *scales3;
  /// basis function counters
  long ibasis0, ibasis1, ibasis2, ibasis3;

  // 'private' iterator fields

  /// shell counters
  long ishell0, ishell1, ishell2, ishell3;
  long ishell3_max;
  /// number of primitives
  long nprim0, nprim1, nprim2, nprim3;
  /// primitive counters
  long iprim0, iprim1, iprim2, iprim3;
  /// primitive offsets
  long oprim0, oprim1, oprim2, oprim3;
};

#endif  // GBASIS_ITER_GB4_H_
