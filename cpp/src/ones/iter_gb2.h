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
 * @file iter_gb2.h
 * @brief Iterators over Gaussian basis functions
 */

#ifndef GBASIS_ITER_GB2_H_
#define GBASIS_ITER_GB2_H_

#include "../gbasis.h"

//! A 2-basis function iterator. Used for 1-electron integrals.
class IterGB2 {
 private:
  // input fields
  const GBasis *gbasis;
  const long *basis_offsets;

 public:
  explicit IterGB2(GBasis *gbasis);
  IterGB2(const IterGB2 &other) = delete;

  /// Increment the shell counters (oprims, ishells) Takes into account 2-fold symmetry.
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

  /**
   * Dot product of the current shell, stored in work, with dm.
   *
   * @param work array of pure basis functions for the current shell.
   * @param dm array which will be dotted with the current shell.
   *
   * @return the dot product of the current shell with the array dm.
   */

  double dot(const double *work, const double *dm);

  // 'public' iterator fields

  /// Shell types
  long shell_type0, shell_type1;
  /// Contraction coefficient
  double con_coeff;
  //@{
  /// Exponents
  double alpha0;
  double alpha1;
  //@}
  //@{
  /// Current gaussian centers
  const double *r0;
  const double *r1;
  //@}
  //@{
  /// Normalization constants
  const double *scales0;
  const double *scales1;
  //@}
  //@{
  /// Basis function counters
  long ibasis0;
  long ibasis1;
  //@}

  // 'private' iterator fields
  //@{
  /// Shell counters
  long ishell0;
  long ishell1;
  //@}
  //@{
  /// number of primitives
  long nprim0;
  long nprim1;
  //@}
  //@{
  // primitive counters
  long iprim0;
  long iprim1;
  //@}
  //@{
  //primitive offsets
  long oprim0;
  long oprim1;
  //@}
};

#endif  // GBASIS_ITER_GB2_H_
