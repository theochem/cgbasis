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
 * @file iter_gb1.h
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
  //! A 1-basis function iterator. Just used for basic structure.
  explicit IterGB1(GBasis *gbasis);
  IterGB1(const IterGB1 &other) = delete;

  //! Increments the shell index
  int inc_shell();

  //! Update internal shell counter and resets center
  void update_shell();

  //! Increment next primitive
  int inc_prim();

  //! Update internal primitive counter
  void update_prim();

  /**
   * Store data into
   * @param work
   * @param output
   * @param dim
   */
  void store(const double *work, double *output, long dim);

  // 'public' iterator fields
  long shell_type0; ///< shell angular momentum
  double con_coeff; ///< contraction coefficient
  double alpha0; ///< shell exponent
  const double *r0;  ///< The current center
  const double *scales0;  ///< Normalization constants
  long ibasis0;  ///< Basis function counters (for output storage)

  // 'private' iterator fields
  long ishell0;  ///< Shell counters
  long nprim0; ///< number of primitives
  long iprim0; ///< primitive counter
  long oprim0;  ///< primitive offset
};

#endif  // GBASIS_ITER_GB1_H_
