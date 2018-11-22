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
 * @file gbasis.h
 * @brief Gaussian basis set classes
 */

#ifndef GBASIS_GBASIS_H_
#define GBASIS_GBASIS_H_

#include "ones/ints2.h"
#include "twos/ints4.h"
#include "grids/fns.h"

/**
 * Normalization constant for cartesian gaussian basis functions
 * @param alpha Exponent
 * @param n Cartesian center
 * @return constant
 */
const double gob_cart_normalization(const double alpha, const long *n);

/**
 * Normalization constant for pure gaussian basis functions
 * @param alpha exponent
 * @param l angular momentum
 * @return constant
 */
const double gob_pure_normalization(const double alpha, const long l);

/**
 * Base class for Gaussian basis integrals
 */
class GBasis {
 private:
  /// Auxiliary arrays that contain convenient derived information.
  long *basis_offsets; ///< offsets of basis functions by center
  long *prim_offsets; ///< offsets for primitives by center
  long *scales_offsets; ///< offsets for normalization constants by center
  long *shell_lookup; ///< shell lookup table

  double *scales;  ///< pre-computed normalization constants.
  long nbasis; ///< number of basis functions
  long nscales; ///< number of normalization constants
  long max_shell_type; ///< maximum angular momentum

 public:
  // Arrays that fully describe the basis set.
  const double *centers; ///< gaussian centers
  const long *shell_map; ///< the center index for each shell
  const long *nprims; ///< number of primitives per shell
  const long *shell_types; ///< shell angular momenta
  const double *alphas; ///< exponents
  const double *con_coeffs; ///< contraction coefficients
  const long ncenter; ///< number of centers
  const long nshell; ///< number of shells
  const long nprim_total; ///< total number of primitives

  double r0[3]; ///< center 0
  double r1[3]; ///< center 1
  double r2[3]; ///< center 2
  double r3[3]; ///< center 3

  /**
   * This class describes basis sets applied to a certain molecular structure.
   *
   * The order of the pure shells is based on the order of real spherical.
   * The functions are sorted from low to high magnetic quantum number,
   * with cosine-like functions before the sine-like functions. The order
   * of functions in a Cartesian shell is alphabetic. Some examples:
   *
   * shell_type = 0, S:
   *  0 -> 1
   * shell_type = 1, P:
   *  0 -> x
   *  1 -> y
   *  2 -> z
   * shell_type = 2, Cartesian D:
   *  0 -> xx
   *  1 -> xy
   *  2 -> xz
   *  3 -> yy
   *  4 -> yz
   *  5 -> zz
   * shell_type = 3, Cartesian F:
   *  0 -> xxx
   *  1 -> xxy
   *  2 -> xxz
   *  3 -> xyy
   *  4 -> xyz
   *  5 -> xzz
   *  6 -> yyy
   *  7 -> yyz
   *  8 -> yzz
   *  9 -> zzz
   * shell_type = -1, not allowed
   * shell_type = -2, pure D:
   *  0 -> zz
   *  1 -> xz
   *  2 -> yz
   *  3 -> xx-yy
   *  4 -> xy
   * shell_type = -3, pure F:
   *  0 -> zzz
   *  1 -> xzz
   *  2 -> yzz
   *  3 -> xxz-yyz
   *  4 -> xyz
   *  5 -> xxx-3xyy
   *  6 -> 3xxy-yyy
   *
   * @param centers centers for the basis functions.
   * @param shell_map the center index for each shell.
   * @param nprims The number of primitives in each shell.
   * @param shell_types contraction types for each shell
   * @param alphas The exponents of the primitives in one shell.
   * @param con_coeffs
   *    The contraction coefficients of the primitives for each
   *     contraction in a contiguous array. The coefficients are ordered
   *     according to the shells. Within each shell, the coefficients are
   *     grouped per exponent.
   * @param ncenter number of centers
   * @param nshell number of shells
   * @param nprim_total total number of primitives
   */
  GBasis(const double *centers, const long *shell_map, const long *nprims,
         const long *shell_types, const double *alphas, const double *con_coeffs,
         const long ncenter, const long nshell, const long nprim_total);
  GBasis(const GBasis &other) = delete;

  virtual ~GBasis();

  virtual const double normalization(const double alpha, const long *n) const = 0;

  void init_scales();

  void shift_center(const double *r, double *shift, double *r_total);

  void compute_two_index(double *output, GB2Integral *integral);

  void compute_four_index(double *output, GB4Integral *integral, double *shift = nullptr);

  void compute_grid_point1(double *output, double *point, GB1GridFn *grid_fn);

  double compute_grid_point2(double *dm, double *point, GB2DMGridFn *grid_fn);

  const long get_nbasis() const { return nbasis; }

  const long get_nscales() const { return nscales; }

  const long get_max_shell_type() const { return max_shell_type; }

  const long *get_basis_offsets() const { return basis_offsets; }

  const long *get_prim_offsets() const { return prim_offsets; }

  const long *get_shell_lookup() const { return shell_lookup; }

  const double *get_scales(long iprim) const { return scales + scales_offsets[iprim]; }
};

class GOBasis : public GBasis {
 public:
  // double compute_normalization(double alpha, long nx, long ny, long nz);
  GOBasis(const double *centers, const long *shell_map, const long *nprims,
          const long *shell_types, const double *alphas, const double *con_coeffs,
          const long ncenter, const long nshell, const long nprim_total);

  const double normalization(const double alpha, const long *n) const;

  /** @brief
          Computes the overlap integrals.

      @param output
          The output array with the integrals.
   */
  void compute_overlap(double *output);

  /** @brief
          Computes the kinetic integrals.

      @param output
          The output array with the integrals.
   */
  void compute_kinetic(double *output);

  /** @brief
          Computes the nuclear attraction integrals.

      @param charges
          The array with values on the nuclear charges.

      @param centers
          The array with location of the nuclear charges.

      @param ncharge
          The number of nuclear charges.

      @param output
          The output array with the integrals.
   */
  void compute_nuclear_attraction(double *charges, double *centers, long ncharge,
                                  double *output);

  /** @brief
          Computes the nuclear attraction integrals.

      @param charges
          The array with values on the nuclear charges.

      @param centers
          The array with location of the nuclear charges.

      @param ncharge
          The number of nuclear charges.

      @param output
          The output array with the integrals.

      @param mu
          The range-separation parameter.
   */
  void compute_erf_attraction(double *charges, double *centers, long ncharge,
                              double *output, double mu);

  /** @brief
          Computes the nuclear attraction integrals.

      @param charges
          The array with values on the nuclear charges.

      @param centers
          The array with location of the nuclear charges.

      @param ncharge
          The number of nuclear charges.

      @param output
          The output array with the integrals.
      @param c
          Coefficient of the gaussian.

      @param alpha
          Exponential parameter of the gaussian.
   */
  void compute_gauss_attraction(double *charges, double *centers, long ncharge,
                                double *output, double c, double alpha);

  /** @brief
          Computes the electron repulsion integrals.

      @param output
          The output array with the integrals.
   */
  void compute_electron_repulsion(double *output);

  /** @brief
          Computes the ERF electron repulsion integrals.

      @param output
          The output array with the integrals.

      @param mu
          The range-separation parameter.
   */
  void compute_erf_repulsion(double *output, double mu);

  /** @brief
          Computes the Gaussian electron repulsion integrals.

      @param output
          The output array with the integrals.

      @param c
          Coefficient of the gaussian.

      @param alpha
          Exponential parameter of the gaussian.
   */
  void compute_gauss_repulsion(double *output, double c, double alpha);

  /** @brief
          Computes the r^alpha electron repulsion integrals.

      @param output
          The output array with the integrals.

      @param alpha
          The power of r in the potential.
   */
  void compute_ralpha_repulsion(double *output, double alpha);

  /** @brief
          Computes the Dirac-Delta electron repulsion integrals.

      @param output
          The output array with the integrals.
   */
  void compute_delta_repulsion(double *output, double *shift = nullptr);

  /** @brief
          Computes the Intracule core integrals.

      @param output
          The output array with the integrals.

      @param point
          The intracular coordinate.
   */
  void compute_intra_density(double *output, double *point);

  /** @brief
          Computes the (multipole) moment integrals.

      @param xyz
          The powers of xyz in the integrals.

      @param center
          The location around which the moment integrals are computed.

      @param output
          The output array with the integrals.
   */
  void compute_multipole_moment(long *xyz, double *center, double *output);

  void compute_grid1_exp(long nfn, double *coeffs, long npoint, double *points,
                         long norb, long *iorbs, double *output);

  /** @brief
          Computes the gradient of the molecular orbital on a grid.

      @param nfn
          The number of functions.

      @param coeffs
          The coefficients for the basisfunction expanion.

      @param npoint
          The number of grid points to be calculated.

      @param points
          The coordinates of grid points to be calculated.

      @param norb
          The number of orbitals to be calculated.

      @param iorbs
          The orbitals to be calculated.

      @param output
          The output array with the integrals.
*/
  void compute_grid1_grad_exp(long nfn, double *coeffs, long npoint,
                              double *points, long norb, long *iorbs, double *output);

  void compute_grid1_dm(double *dm, long npoint, double *points,
                        GB1DMGridFn *grid_fn, double *output,
                        double epsilon, double *dmmaxrow);

  void compute_grid2_dm(double *dm, long npoint, double *points, double *output);

  void compute_grid1_fock(long npoint, double *points, double *weights,
                          long pot_stride, double *pots,
                          GB1DMGridFn *grid_fn, double *output);
};

#endif  // GBASIS_GBASIS_H_
