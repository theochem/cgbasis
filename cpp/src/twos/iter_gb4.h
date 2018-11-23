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
 * @file iter_gb4.h
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

  long shell_type0; ///< shell types
  long shell_type1; ///< shell types
  long shell_type2; ///< shell types
  long shell_type3; ///< shell types
  /// Contraction coefficient
  double con_coeff;

  double alpha0; ///< Exponents
  double alpha1; ///< Exponents
  double alpha2; ///< Exponents
  double alpha3; ///< Exponents

  const double *r0; ///< current gaussian centers
  const double *r1; ///< current gaussian centers
  const double *r2; ///< current gaussian centers
  const double *r3; ///< current gaussian centers

  const double *scales0; ///< normalization constants
  const double *scales1; ///< normalization constants
  const double *scales2; ///< normalization constants
  const double *scales3; ///< normalization constants

  long ibasis0; ///< basis function counters
  long ibasis1; ///< basis function counters
  long ibasis2; ///< basis function counters
  long ibasis3; ///< basis function counters

  // 'private' iterator fields

  long ishell0; ///< shell counters
  long ishell1; ///< shell counters
  long ishell2; ///< shell counters
  long ishell3; ///< shell counters

  long ishell3_max; ///< shell

  long nprim0; ///< number of primitives
  long nprim1; ///< number of primitives
  long nprim2; ///< number of primitives
  long nprim3; ///< number of primitives

  long iprim0; ///< primitive counters
  long iprim1; ///< primitive counters
  long iprim2; ///< primitive counters
  long iprim3; ///< primitive counters

  long oprim0; ///< primitive offsets
  long oprim1; ///< primitive offsets
  long oprim2; ///< primitive offsets
  long oprim3; ///< primitive offsets
};

#endif  // GBASIS_ITER_GB4_H_
