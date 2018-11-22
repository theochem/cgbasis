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
 * @file calc.h
 * @brief Base class for any integral/evaluation of Gaussian functions
 */

#ifndef GBASIS_CALC_H_
#define GBASIS_CALC_H_

/** @brief
      Base class for anything that computes stuff by iterating over Gaussian
      basis functions.

    */
class GBCalculator {
 protected:
  long nwork; ///< number of elements in work array
  long max_shell_type; ///< maximum shell angular momentum
  long max_nbasis; ///< maximum number of basis functions
  double *work_pure; ///< contiguous pure work arrays sufficiently large for max_shell_type
  double *work_cart;  ///< contiguous cartesian work arrays sufficiently large for max_shell_type
  void swap_work(); ///< write work_pure into work_cart

 public:
  /** @brief
        Construct a GBCalculator object.

      This also allocates work arrays for manipulating and storing intermediate
      results. The size of these work arrays is dim_work*max_nbasis**basis_work.

      @param max_shell_type
        The maximum shell type in the basis set. This is used to allocate
        sufficiently large working arrays.

      @param dim_work
        Prefactor for the size of the work arrays.

      @param basis_work
        The work array size is multiplied by max_nbasis**basis_work.
    */
  GBCalculator(long max_shell_type, long dim_work, int basis_work);
  GBCalculator(const GBCalculator& other) = delete;

  virtual ~GBCalculator();
  /** @brief
    Number of elements in the work array.
    */
  const long get_nwork() const { return nwork; }
/** @brief
    Maximum shell angular momentum
    */
  const long get_max_shell_type() const { return max_shell_type; }
/** @brief
    Maximum number of basis functions
    */
  const long get_max_nbasis() const { return max_nbasis; }
/** @brief
    The cartesian work array
    */
  const double *get_work() const { return work_cart; }
};

#endif  // GBASIS_CALC_H_
