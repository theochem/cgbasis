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
 * @file ints4.h
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

#ifndef GBASIS_INTS4_H_
#define GBASIS_INTS4_H_

#include "libint2.h"
#include "../calc.h"

//! Base class for four-center integrals.
class GB4Integral : public GBCalculator {
 public:
  /** @brief
          Initialize a GB4Integral object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.
    */
  explicit GB4Integral(long max_shell_type)
    : GBCalculator(max_shell_type, 1, 4), shell_type0(0), shell_type1(0),
      shell_type2(0), shell_type3(0), r0(NULL), r1(NULL), r2(NULL), r3(NULL) {}

  /** @brief
          Set internal parameters for a new group of four contractions.

      @param shell_type0
          Angular momentum index for contraction 0.

      @param shell_type1
          Angular momentum index for contraction 1.

      @param shell_type2
          Angular momentum index for contraction 2.

      @param shell_type3
          Angular momentum index for contraction 3.

      @param r0
          Cartesian coordinates of center 0.

      @param r1
          Cartesian coordinates of center 1.

      @param r2
          Cartesian coordinates of center 2.

      @param r3
          Cartesian coordinates of center 3.
    */
  virtual void reset(long shell_type0, long shell_type1, long shell_type2,
                     long shell_type3, const double *r0, const double *r1,
                     const double *r2, const double *r3);

  /** @brief
          Add results for a combination of Cartesian primitive shells to the work array.

      @param coeff
          Product of the contraction coefficients of the four primitives.

      @param alpha0
          The exponent of primitive shell 0.

      @param alpha1
          The exponent of primitive shell 1.

      @param alpha2
          The exponent of primitive shell 2.

      @param alpha3
          The exponent of primitive shell 3.

      @param scales0
          The normalization prefactors for basis functions in primitive shell 0

      @param scales1
          The normalization prefactors for basis functions in primitive shell 1

      @param scales2
          The normalization prefactors for basis functions in primitive shell 2

      @param scales3
          The normalization prefactors for basis functions in primitive shell 3
    */
  virtual void add(double coeff, double alpha0, double alpha1, double alpha2,
                   double alpha3, const double *scales0, const double *scales1,
                   const double *scales2, const double *scales3) = 0;

  /// Transform the results in the work array from Cartesian to pure functions where needed.
  void cart_to_pure();

  /// Shell type of contraction 0
  const long get_shell_type0() const { return shell_type0; }
  /// Shell type of contraction 1
  const long get_shell_type1() const { return shell_type1; }
  /// Shell type of contraction 2
  const long get_shell_type2() const { return shell_type2; }
  /// Shell type of contraction 3
  const long get_shell_type3() const { return shell_type3; }

 protected:
  long shell_type0;  ///< Shell type of contraction 0
  long shell_type1;  ///< Shell type of contraction 1
  long shell_type2;  ///< Shell type of contraction 2
  long shell_type3;  ///< Shell type of contraction 3
  const double *r0;  ///< Center of contraction 0
  const double *r1;  ///< Center of contraction 1
  const double *r2;  ///< Center of contraction 2
  const double *r3;  ///< Center of contraction 3
};

//! Arguments associated with one primitive shell in LibInt conventions.
typedef struct {
  unsigned int am;  ///< Shell type in LibInt conventions.
  const double *r;  ///< Center of a primitive shell.
  double alpha;     ///< Exponent of a primitive shell.
} libint_arg_t;

//! Base class for four-center integrals that use LibInt.
class GB4IntegralLibInt : public GB4Integral {
 public:
  /** @brief
          Initialize a GB4IntegralLibInt object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.
    */
  explicit GB4IntegralLibInt(long max_shell_type);

  ~GB4IntegralLibInt();

  /** @brief
          Set internal parameters for a new group of four contractions.

      See base class for details.
    */
  virtual void reset(long shell_type0, long shell_type1, long shell_type2, long shell_type3,
                     const double *r0, const double *r1, const double *r2, const double *r3);

  /** @brief
          Add results for a combination of Cartesian primitive shells to the work array.

      See base class for details.
    */
  virtual void add(double coeff, double alpha0, double alpha1, double alpha2, double alpha3,
                   const double *scales0, const double *scales1, const double *scales2,
                   const double *scales3);

  /** @brief
          Evaluate the Laplace transform of the the potential.

      For theoretical details and the precise definition of the Laplace transform, we
      refer to the following paper:

      Ahlrichs, R. A simple algebraic derivation of the Obara-Saika scheme for general
      two-electron interaction potentials. Phys. Chem. Chem. Phys. 8, 3072–3077 (2006).
      10.1039/B605188J

      For the general definition of this transform, see Eq. (8) in the reference above.
      Section 5 contains solutions of the Laplace transform for several popular cases.

      @param prefac
          Prefactor with which all results in the output array are multiplied.

      @param rho
          See Eq. (3) in Ahlrichs' paper.

      @param t
          Rescaled distance between the two centers obtained from the application of the
          Gaussian product theorem. See Eq. (5) in Ahlrichs' paper.

      @param mmax
          Maximum derivative of the Laplace transform to be considered.

      @param output
          Output array. The size must be at least mmax + 1.
   */
  virtual void laplace_of_potential(double prefac, double rho, double t, long mmax,
                                    double *output) = 0;

 private:
  Libint_eri_t erieval;         // !< LibInt runtime object.
  libint_arg_t libint_args[4];  // !< Arguments (shell info) for libint.
  long order[4];                // !< Re-ordering of shells for compatibility with LibInt.
  double ab[3];                 // !< Relative vector from shell 2 to 0 (LibInt order).
  double cd[3];                 // !< Relative vector from shell 3 to 1 (LibInt order).
  double ab2;                   // !< Norm squared of ab.
  double cd2;                   // !< Norm squared of cd.
};

/** @brief
        Electron repulsion four-center integrals.

    The potential is 1/r.
  */
class GB4ElectronRepulsionIntegralLibInt : public GB4IntegralLibInt {
 public:
  /** @brief
          Initialize a GB4ElectronRepulsionIntegralLibInt object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.
    */
  explicit GB4ElectronRepulsionIntegralLibInt(long max_shell_type)
      : GB4IntegralLibInt(max_shell_type) {}

  /** @brief
          Evaluate the Laplace transform of the ordinary Coulomb potential.

      See Eq. (39) in Ahlrichs' paper. This is basically a rescaled Boys function.

      See base class for more details.
    */
  virtual void laplace_of_potential(double prefac, double rho, double t, long mmax,
                                    double *output);
};

/** @brief
        Short-range electron repulsion four-center integrals.

    The potential is erf(mu*r)/r.
  */
class GB4ErfIntegralLibInt : public GB4IntegralLibInt {
 public:
  /** @brief
          Initialize a GB4ErfIntegralLibInt object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.

      @param mu
          The range-separation parameter
    */
  GB4ErfIntegralLibInt(long max_shell_type, double mu)
      : GB4IntegralLibInt(max_shell_type), mu(mu) {}

  /** @brief
          Evaluate the Laplace transform of the long-range Coulomb potential.
          (The short-range part is damped away using an error function.) See (52) in
          Ahlrichs' paper.

      See base class for more details.
    */
  virtual void laplace_of_potential(double prefac, double rho, double t, long mmax,
                                    double *output);

  /// The range-separation parameter.
  const double get_mu() const { return mu; }

 private:
  double mu;  // !< The range-separation parameter.
};

/** @brief
        Gaussian electron repulsion four-center integrals.

    The potential is c exp(-alpha r^2).
  */
class GB4GaussIntegralLibInt : public GB4IntegralLibInt {
 public:
  /** @brief
          Initialize a GB4GaussIntegralLibInt object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.

      @param c
          Coefficient of the gaussian.

      @param alpha
          Exponential parameter of the gaussian.
    */
  GB4GaussIntegralLibInt(long max_shell_type, double c, double alpha)
      : GB4IntegralLibInt(max_shell_type), c(c), alpha(alpha) {}

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
  virtual void laplace_of_potential(double prefac, double rho, double t, long mmax,
                                    double *output);

  /// Coefficient of the gaussian.
  const double get_c() const { return c; }
  /// Exponential parameter of the gaussian.
  const double get_alpha() const { return alpha; }

 private:
  double c;  // !< Coefficient of the gaussian.
  double alpha;  // !< Exponential parameter of the gaussian.
};

/** @brief
        Gaussian electron repulsion four-center integrals.

    The potential is r^alpha.
  */
class GB4RAlphaIntegralLibInt : public GB4IntegralLibInt {
 public:
  /** @brief
          Initialize a GB4RAlphaIntegralLibInt object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.

      @param alpha
          The power of r in the potential.
    */
  GB4RAlphaIntegralLibInt(long max_shell_type, double alpha)
      : GB4IntegralLibInt(max_shell_type), alpha(alpha) {}

  /** @brief
          Evaluate the Laplace transform of the r^alpha potential. See Eq. (49) in
          Ahlrichs' paper.

      See base class for more details.
    */
  virtual void laplace_of_potential(double prefac, double rho, double t, long mmax,
                                    double *output);

  /// The power of r.
  const double get_alpha() const { return alpha; }

 private:
  double alpha;  // !< The power of r.
};

/** @brief
        Delta repulsion integrals.

    The potential is \delta(r).
  */
class GB4DeltaIntegralLibInt : public GB4IntegralLibInt {
 public:
  /** @brief
          Initialize a GB4IntraDensIntegralLibInt object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.
    */
  explicit GB4DeltaIntegralLibInt(long max_shell_type)
      : GB4IntegralLibInt(max_shell_type) {}

  /** @brief
          Evaluate the Laplace transform of the ordinary Coulomb potential.

      See Eq. (39) in Ahlrichs' paper. This is basically a rescaled Boys function.

      See base class for more details.
    */
  virtual void laplace_of_potential(double prefac, double rho, double t, long mmax,
                                    double* output);
};


//! Base class for four-center density integrals that use LibInt.
class GB4DIntegralLibInt : public GB4Integral {
 public:
  /** @brief
          Initialize a GB4DensIntegralLibInt object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.
    */
  explicit GB4DIntegralLibInt(long max_shell_type);
  ~GB4DIntegralLibInt();

  /** @brief
          Set internal parameters for a new group of four contractions.

      See base class for details.
    */
  virtual void reset(long shell_type0, long shell_type1, long shell_type2, long shell_type3,
                     const double* r0, const double* r1, const double* r2, const double* r3);
  /** @brief
          Add results for a combination of Cartesian primitive shells to the work array.

      See base class for details.
    */
  virtual void add(double coeff, double alpha0, double alpha1, double alpha2, double alpha3,
                   const double* scales0, const double* scales1, const double* scales2,
                   const double* scales3);

  /** @brief
          Evaluate the Laplace transform of the the potential.

      For theoretical details and the precise definition of the Laplace transform, we
      refer to the following paper:

      Ahlrichs, R. A simple algebraic derivation of the Obara-Saika scheme for general
      two-electron interaction potentials. Phys. Chem. Chem. Phys. 8, 3072–3077 (2006).
      10.1039/B605188J

      For the general definition of this transform, see Eq. (8) in the reference above.
      Section 5 contains solutions of the Laplace transform for several popular cases.

      @param prefac
          Prefactor with which all results in the output array are multiplied.

      @param rho
          See Eq. (3) in Ahlrichs' paper.

      @param t
          Rescaled distance between the two centers obtained from the application of the
          Gaussian product theorem. See Eq. (5) in Ahlrichs' paper.

      @param mmax
          Maximum derivative of the Laplace transform to be considered.

      @param output
          Output array. The size must be at least mmax + 1.
   */
  virtual void laplace_of_potential(double prefac, double rho, double t, double* p,
                                    double* q, long mmax, double* output) = 0;

 private:
  Libint_eri_t erieval;         //!< LibInt runtime object.
  libint_arg_t libint_args[4];  //!< Arguments (shell info) for libint.
  long order[4];                //!< Re-ordering of shells for compatibility with LibInt.
  double ab[3];                 //!< Relative vector from shell 2 to 0 (LibInt order).
  double cd[3];                 //!< Relative vector from shell 3 to 1 (LibInt order).
  double ab2;                   //!< Norm squared of ab.
  double cd2;                   //!< Norm squared of cd.
};

/** @brief
        Intracular (four-center) density integrals at coordinate point.

    The potential is \delta(r - point).
  */
class GB4IntraDensIntegralLibInt : public GB4DIntegralLibInt {
 public:
  /** @brief
          Initialize a GB4IntraDensIntegralLibInt object.

      @param max_shell_type
          Highest angular momentum index to be expected in the reset method.
    */
  GB4IntraDensIntegralLibInt(long max_shell_type, double* point)
      : GB4DIntegralLibInt(max_shell_type), point(point) {}

  /** @brief
          Evaluate the Laplace transform of the ordinary Coulomb potential.

      See Eq. (39) in Ahlrichs' paper. This is basically a rescaled Boys function.

      See base class for more details.
    */
  virtual void laplace_of_potential(double prefac, double rho, double t, double* p,
                                    double* q, long mmax, double* output);

 private:
  double* point;    //!< Array with point.
};

#endif  // GBASIS_INTS4_H_
