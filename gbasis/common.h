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
 * @file common.h
 * @brief Auxiliary functions
 */

#ifndef GBASIS_COMMON_H_
#define GBASIS_COMMON_H_

/// Maximum angular momentum this library has been configured for.
#define MAX_SHELL_TYPE 7
/// Maximum size of work arrays in Cartesian polynomials
#define MAX_NCART_CUMUL ((MAX_SHELL_TYPE+1)*(MAX_SHELL_TYPE+2)*(MAX_SHELL_TYPE+3))/6
/// Maximum size of derivative work arrays in Cartesian polynomials
#define MAX_NCART_CUMUL_D ((MAX_SHELL_TYPE+2)*(MAX_SHELL_TYPE+3)*(MAX_SHELL_TYPE+4))/6
/// Maximum size of 2nd derivative work arrays in Cartesian polynomials
#define MAX_NCART_CUMUL_DD ((MAX_SHELL_TYPE+3)*(MAX_SHELL_TYPE+4)*(MAX_SHELL_TYPE+5))/6

// Simple math stuff
/// Factorial
long fac(long n);

/**
 * Factorial of every other number
 * \f[ n(n-2)(n-4)(n-8)...1 \f]
 * @param n Upper value
 */
long fac2(long n);

/// Binomial coefficient
long binom(long n, long m);

/**
 * Get number of basis functions in shell
 *
 * @param shell_type 0=S (cart), 1=P (cart), 2=D (cart), -1=S (pure), -2=D (pure), etc.
 */
long get_shell_nbasis(long shell_type);


/**
 * Return the largest shell that Libint (and subsequently, GBasis) has been configured for.
 * @return 0=S, 1=P, 2=D, and so on.
 */
long get_max_shell_type();

const double dist_sq(const double *r0, const double *r1);

// Auxiliary functions for Gaussian integrals
/**
 * The sum of two Gaussian centres (via the gaussian product theorem).
 *
 * \f[ \frac{\alpha_{1}R_{0}+\alpha_{2}R_{1}}{\alpha_{1}+\alpha_{2}} \f]
 *
 * See eqn 2.3 and preceeding line in Taketa et al. (1966)
 *
 * @param alpha0 Exponent on centre 0
 * @param r0 Centre 0
 * @param alpha1 Exponent on centre 1
 * @param r1 Centre 1
 * @param gamma_inv \f$ \frac{1}{\alpha_{1}+\alpha_{2}} \f$
 * @param gpt_center Output variable. The resultant gaussian centre.
 */
void compute_gpt_center(double alpha0, const double *r0, double alpha1, const double *r1,
                        double gamma_inv,
                        double *gpt_center);

/**
 * Gaussian product theorem coefficient
 *
 * \f[ f_{k}(n_{0},n_{1},pa,pb) \f]
 *
 * See eqn 2.4 in Taketa et al. (1966)
 *
 * @param k
 * @param n0 angular momentum quantum number of primitive 1D shell 0
 * @param n1 angular momentum quantum number of primitive 1D shell 1
 * @param pa gaussian centre of shell 0
 * @param pb gaussian centre of shell 1
 */
double gpt_coeff(long k, long n0, long n1, double pa, double pb);

/**
 * Gaussian basis overlap integral (1D)
 *
 * \f[ \sum_{k}^{(n0+n1)/2}f_{2k}(n_{0},n_{1},pa,pb)\frac{(2k-1)!!}{(2\gamma)^{k}} \f]
 *
 * See eqn 2.12 in Taketa et al. (1966)
 *
 * @param n0 angular momentum quantum number of primitive 1D shell 0
 * @param n1 angular momentum quantum number of primitive 1D shell 1
 * @param pa gaussian centre 0
 * @param pb gaussian centre 1
 * @param gamma_inv \f$ \frac{1}{\alpha_1 + \alpha_2}\f$
 */
double gb_overlap_int1d(long n0, long n1, double pa, double pb, double gamma_inv);

void nuclear_attraction_helper(double *work_g, long n0, long n1, double pa, double pb, double pc,
                               double gamma_inv);

// Auxiliary functions for r^alpha integrals

//! \f[ cit=t^{m}\frac{2^{2i}}{(2i+1)!} \f]
double cit(int i, double t, int m);

//! \f[ j(j-1)(j-2)(j-3)...(j-n) \f]
long jfac(int j, int n);

/** @brief Evaluate the taylor series for r^alpha integrals

    \f[ \sum_{i}t^{i}\Gamma(i+\frac{\alpha+3}{2})2^{\frac{2i}{(2i+1)!}} \f]

    @param n
        Angular moment (the order of the derivative of the basic integral Gn in Alhrichs
        Phys. Chem. Chem. Phys., 8, 3072 (2006)). The maximum value implemented is n=10.

    @param alpha
        The power of r in the potential.

    @param t
        \f$ \rho|p-q|^{2} \f$

    @param prefac
        \f$ \frac{e^{-t}}{\rho^{\frac{3}{2}}} \f$ - This term helps the Taylor series to converge when t is a large
        number, the factor \f$ \frac{1}{2}\sqrt{\rho^{\alpha}} \f$ was "replaced" and multiplied outside, at
        the end, in the laplace_of_potential function.
*/
double dtaylor(int n, double alpha, double t, double prefac);

#endif  // GBASIS_COMMON_H_
