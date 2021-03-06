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
 * @file ints2.h
 * @brief Evaluation of integrals of Gaussian basis functions
 *
    The implementation of the two-index operators in this module are based on
    the paper "Gaussian-Expansion Methods for Molecular Integrals", H. Taketa,
    S. Huzinga, K. O-ohata, Journal of the Physical Society of Japan, vol. 21,
    p. 2313, y. 1966. Be aware that there are some misprints in the paper:

    - eq. 2.18: CP_x should be p_x, i.e. change of sign.
    - list of expressions at the bottom of p. 2319:
        case l_1+l_2=2, G_1 = -f_1*p - f_2/(2*gamma)

    The implementation of the (multipole) moment integrals in this module are
    based on the paper "Efficient recursive computation of molecular integrals
    over Cartesian Gaussian functions", S. Obara, A. Saika, Journal of Chemical
    Physics, vol. 84, p. 3963, y. 1986.
 *
 */

#ifndef GBASIS_INTS2_H_
#define GBASIS_INTS2_H_

#include <cstddef>
#include "../calc.h"
#include "iter_pow2.h"

/**
 * Integral of two basis functions. Commonly used for 1-electron integrals.
 */
class GB2Integral : public GBCalculator {
 protected:
  long shell_type0; ///< shell type 0
  long shell_type1; ///< shell type 1
  const double *r0; ///< gaussian center 0
  const double *r1; ///< Gaussian center 1
  IterPow2 i2p; ///< Cartesian basis function iterator
 public:
  /**
   * Construct a GB2Integral object.
   *
   * @param max_shell_type The maximum shell type in the basis set. This is used to allocate
   *    sufficiently large working arrays.
   */
  explicit GB2Integral(long max_shell_type)
      : GBCalculator(max_shell_type, 1, 2), shell_type0(0), shell_type1(0),
        r0(NULL), r1(NULL), i2p() {}
  /**
   * Re-initialize on a new shell.
   *
   * @param shell_type0 Shell 0 type
   * @param shell_type1 Shell 1 type
   * @param r0 Gaussian Centre 0
   * @param r1 Gaussian Centre 1
   */
  void reset(long shell_type0, long shell_type1, const double *r0, const double *r1);

  /**
   Add results for a combination of Cartesian primitive shells to the work array.

   See eqn 2.12 of Taketa et al. (1966)

   @param coeff
       Product of the contraction coefficients of the two primitives.

   @param alpha0
       The exponent of primitive shell 0.

   @param alpha1
       The exponent of primitive shell 1.

   @param scales0
      The normalization prefactors for basis functions in primitive shell 0

   @param scales1
      The normalization prefactors for basis functions in primitive shell 1
   */
  virtual void add(double coeff, double alpha0, double alpha1, const double *scales0,
                   const double *scales1) = 0;

  /**
   * Transform the work arrays from cartesian into pure coordinates. N.B. The results are store back into work_cart!
   */
  void cart_to_pure();

  /**
   * Angular momentum of shell 0
   * @return
   */
  const long get_shell_type0() const { return shell_type0; }

  /**
   * Angular momentum of shell 1
   * @return
   */
  const long get_shell_type1() const { return shell_type1; }
};

/** @brief
 Compute the overlap integrals in a Gaussian orbital basis.
 */
class GB2OverlapIntegral : public GB2Integral {
 public:
  /**
   * Compute the overlap integrals in a Gaussian orbital basis.
   * @param max_shell_type maximum shell type of integral
   */
  explicit GB2OverlapIntegral(long max_shell_type) : GB2Integral(max_shell_type) {}

  virtual void
  add(double coeff, double alpha0, double alpha1, const double *scales0, const double *scales1);
};

/** @brief
 Compute the kinetic integrals in a Gaussian orbital basis.
 */
class GB2KineticIntegral : public GB2Integral {
 public:
  /**
   * Compute the kinetic integrals in a Gaussian orbital basis.
   * @param max_shell_type maximum shell type of integral
   */
  explicit GB2KineticIntegral(long max_shell_type) : GB2Integral(max_shell_type) {}

  virtual void
  add(double coeff, double alpha0, double alpha1, const double *scales0, const double *scales1);
};

/** @brief
 Compute the nuclear attraction integrals in a Gaussian orbital basis.
 */
class GB2AttractionIntegral : public GB2Integral {
 private:
  double *charges;    // !< Array with values of the nuclear charges.
  double *centers;    // !< The centers where the charges are located.
  long ncharge;       // !< Number of nuclear charges.

  double *work_g0;    // !< Temporary array to store intermediate results.
  double *work_g1;    // !< Temporary array to store intermediate results.
  double *work_g2;    // !< Temporary array to store intermediate results.
  double *work_boys;  // !< Temporary array to store the laplace of the interaction potential.

 public:
  /** @brief
        Initialize a GB2AttractionIntegral object.

   @param max_shell_type
       Highest angular momentum index to be expected in the reset method.

   */
  GB2AttractionIntegral(long max_shell_type, double *charges, double *centers, long ncharge);

  ~GB2AttractionIntegral();

  /** @brief
    Add results for a combination of Cartesian primitive shells to the work array.

   @param coeff
       Product of the contraction coefficients of the two primitives.

   @param alpha0
       The exponent of primitive shell 0.

   @param alpha1
       The exponent of primitive shell 1.

   @param scales0
      The normalization prefactors for basis functions in primitive shell 0

   @param scales1
      The normalization prefactors for basis functions in primitive shell 1
  */
  virtual void add(double coeff, double alpha0, double alpha1, const double *scales0,
                   const double *scales1);

  /** @brief
    Evaluate the Laplace transform of the the potential applied to nuclear attraction terms.

   For theoretical details and the precise definition of the Laplace transform, we
   refer to the following paper:

   Ahlrichs, R. A simple algebraic derivation of the Obara-Saika scheme for general
   two-electron interaction potentials. Phys. Chem. Chem. Phys. 8, 3072–3077 (2006).
   10.1039/B605188J

   For the general definition of this transform, see Eq. (8) in the reference above.
   Section 5 contains solutions of the Laplace transform for several popular cases.

   @param gamma
       Sum of the exponents of the two gaussian functions involved in the integral.
       Similar to  the first term in Eq. (3) in Ahlrichs' paper.

   @param arg
       Rescaled distance between the two centers obtained from the application of the
       Gaussian product theorem. Equivalent to Eq. (5) in Ahlrichs' paper.

   @param mmax
       Maximum derivative of the Laplace transform to be considered.

   @param output
       Output array. The size must be at least mmax + 1.
   */
  virtual void laplace_of_potential(double gamma, double arg, long mmax, double *output) = 0;
};

/** @brief
        Nuclear Electron Attraction two-center integrals.

    The potential is 1/r.
  */
class GB2NuclearAttractionIntegral : public GB2AttractionIntegral {
 public:
  /** @brief
          Initialize a GB2NuclearAttractionIntegral object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.

      @param charges
          Array with values of the charges.

      @param centers
          The centers [[C1_x, C1_y, C1_z],[...],] around which the moment integrals are computed.

      @param ncharge
          Number of nuclear charges.
   */
  explicit GB2NuclearAttractionIntegral(long max_shell_type, double *charges,
                                        double *centers, long ncharge)
      : GB2AttractionIntegral(max_shell_type, charges, centers, ncharge) {}

  /** @brief
          Evaluate the Laplace transform of the ordinary Coulomb potential.

      See Eq. (39) in Ahlrichs' paper. This is basically a rescaled Boys function.

      See base class for more details.
    */
  virtual void laplace_of_potential(double gamma, double arg, long mmax, double *output);
};

/** @brief
        Short-range electron repulsion four-center integrals.

    The potential is erf(mu*r)/r.
  */
class GB2ErfAttractionIntegral : public GB2AttractionIntegral {
 public:
  /** @brief
          Initialize a GB2ErfAttractionIntegral object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.

      @param charges
          Array with values of the charges.

      @param centers
          The centers [[C1_x, C1_y, C1_z],[...],] around which the moment integrals are computed.

      @param ncharge
          Number of nuclear charges.

      @param mu
          The range-separation parameter.
    */
  GB2ErfAttractionIntegral(long max_shell_type, double *charges, double *centers,
                           long ncharge, double mu)
      : GB2AttractionIntegral(max_shell_type, charges, centers, ncharge), mu(mu) {}

  /** @brief
          Evaluate the Laplace transform of the long-range Coulomb potential.
          (The short-range part is damped away using an error function.) See (52) in
          Ahlrichs' paper.

      See base class for more details.
    */
  virtual void laplace_of_potential(double gamma, double arg, long mmax, double *output);

  /**
   * The range-separation parameter.
   * @return
   */
  const double get_mu() const { return mu; }

 private:
  double mu;  // !< The range-separation parameter.
};

/** @brief
        Gaussian nuclear electron attraction two-center integrals.

    The potential is c exp(-alpha r^2).
  */
class GB2GaussAttractionIntegral : public GB2AttractionIntegral {
 public:
  /** @brief
          Initialize a GB2GaussAttractionIntegral object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.

      @param charges
          Array with values of the charges.

      @param centers
          The centers [[C1_x, C1_y, C1_z],[...],] around which the moment integrals are computed.

      @param ncharge
          Number of nuclear charges.

      @param c
          Coefficient of the gaussian.

      @param alpha
          Exponential parameter of the gaussian.
    */
  GB2GaussAttractionIntegral(long max_shell_type, double *charges, double *centers, long ncharge,
                             double c, double alpha)
      : GB2AttractionIntegral(max_shell_type, charges, centers, ncharge), c(c),
        alpha(alpha) {}

  /** @brief
          Evaluate the Laplace transform of the Gaussian potential.

          See Ahlrichs' paper for details. This type of potential is used in the papers
          of P.M.W Gill et al. and J. Toulouse et al.:

          Gill, P. M. W., & Adamson, R. D. (1996). A family of attenuated Coulomb
          operators. Chem. Phys. Lett., 261(1-2), 105–110.
          http://doi.org/10.1016/0009-2614(96)00931-1

          Toulouse, J., Colonna, F., & Savin, A. (2004). Long-range-short-range separation
          of the electron-electron interaction in density-functional theory. Phys. Rev. A,
          70, 62505. http://doi.org/10.1103/PhysRevA.70.062505

      See base class for more details.
    */
  virtual void laplace_of_potential(double gamma, double arg, long mmax, double *output);

  /**
   * Coefficient of the gaussian.
   * @return
   */
  const double get_c() const { return c; }
  /**
   * Exponential parameter of the gaussian.
   * @return
   */
  const double get_alpha() const { return alpha; }

 private:
  double c;  // !< Coefficient of the gaussian.
  double alpha;  // !< Exponential parameter of the gaussian.
};

/** @brief
        Compute the (multipole) moment integrals in a Gaussian orbital basis.
        < gto_a | (x - C_x)^l (y - C_y)^m (z - C_z)^n | gto_b >.
 */
class GB2MomentIntegral : public GB2Integral {
 private:
  long *xyz;          // !< Powers for x, y and z of the multipole moment.
  double *center;     // !< The origin w.r.t. to which the multipole moment is computed.

 public:
  /** @brief
          Initialize Moment integral calculator

      @param max_shell_type
          The highest angular momentum index suported

      @param xyz
          The powers of x,y,z in the integrals (l, m, n).

      @param center
          The center [C_x, C_y, C_z] around which the moment integrals arecomputed
  */
  GB2MomentIntegral(long max_shell_type, long *xyz, double *center);

  /** @brief
          Add integrals for a pair of primite shells to the current contraction.

      @param coeff
          The contraction coefficient for the current primitive.

      @param alpha0
          The exponent of the primitive shell 0.

      @param alpha1
          The exponent of the primitive shell 1.

      @param scales0
          The normalization constants for the basis functions in primitive shell 0.

      @param scales1
          The normalization constants for the basis functions in primitive shell 1.
    */
  virtual void add(double coeff, double alpha0, double alpha1,
                   const double *scales0, const double *scales1);
};

#endif  // GBASIS_INTS2_H_
