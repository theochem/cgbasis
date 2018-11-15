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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include "../twos/boys.h"
#include "../cartpure.h"
#include "ints2.h"

using std::abs;




/*

   GB2Integral

*/


void
GB2Integral::reset(long _shell_type0, long _shell_type1, const double *_r0, const double *_r1) {
  if ((_shell_type0 < -max_shell_type) || (_shell_type0 > max_shell_type)) {
    throw std::domain_error("shell_type0 out of range.");
  }
  if ((_shell_type1 < -max_shell_type) || (_shell_type1 > max_shell_type)) {
    throw std::domain_error("shell_type1 out of range.");
  }
  shell_type0 = _shell_type0;
  shell_type1 = _shell_type1;
  r0 = _r0;
  r1 = _r1;
  // We make use of the fact that a floating point zero consists of
  // consecutive zero bytes.
  memset(work_cart, 0, nwork * sizeof(double));
  memset(work_pure, 0, nwork * sizeof(double));
}

void GB2Integral::cart_to_pure() {
  /*
     The initial results are always stored in work_cart. The projection
     routine always outputs its result in work_pure. Once that is done,
     the pointers to both blocks are swapped such that the final result is
     always back in work_cart.
  */

  // Project along index 0 (rows)
  if (shell_type0 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type0,
                     1,  // anterior
                     get_shell_nbasis(abs(shell_type1)));  // posterior
    swap_work();
  }

  // Project along index 1 (cols)
  if (shell_type1 < -1) {
    cart_to_pure_low(work_cart, work_pure, -shell_type1,
                     get_shell_nbasis(shell_type0),  // anterior
                     1);  // posterior
    swap_work();
  }
}

/*

   GB2OverlapIntegral

*/


void GB2OverlapIntegral::add(double coeff, double alpha0, double alpha1, const double *scales0,
                             const double *scales1) {
  double pre, gamma_inv;
  double gpt_center[3];

  gamma_inv = 1.0 / (alpha0 + alpha1);
  pre = coeff * exp(-alpha0 * alpha1 * gamma_inv * dist_sq(r0, r1));
  compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
  i2p.reset(abs(shell_type0), abs(shell_type1));
  do {
    work_cart[i2p.offset] += pre * (
        gb_overlap_int1d(i2p.n0[0], i2p.n1[0], gpt_center[0] - r0[0], gpt_center[0] - r1[0],
                         gamma_inv) *
            gb_overlap_int1d(i2p.n0[1], i2p.n1[1], gpt_center[1] - r0[1], gpt_center[1] - r1[1],
                             gamma_inv) *
            gb_overlap_int1d(i2p.n0[2], i2p.n1[2], gpt_center[2] - r0[2], gpt_center[2] - r1[2],
                             gamma_inv) *
            scales0[i2p.ibasis0] * scales1[i2p.ibasis1]);
  } while (i2p.inc());
}

/*

   GB2KineticIntegral

*/


double kinetic_helper(double alpha0, double alpha1, long n0, long n1, double pa, double pb,
                      double gamma_inv) {
  double poly = 0;

  if (n0 > 0) {
    if (n1 > 0) {
      // <-1|-1>
      poly += 0.5 * n0 * n1 * gb_overlap_int1d(n0 - 1, n1 - 1, pa, pb, gamma_inv);
    }
    // <-1|+1>
    poly -= alpha1 * n0 * gb_overlap_int1d(n0 - 1, n1 + 1, pa, pb, gamma_inv);
  }
  if (n1 > 0) {
    // <+1|-1>
    poly -= alpha0 * n1 * gb_overlap_int1d(n0 + 1, n1 - 1, pa, pb, gamma_inv);
  }
  // <+1|+1>
  poly += 2.0 * alpha0 * alpha1 * gb_overlap_int1d(n0 + 1, n1 + 1, pa, pb, gamma_inv);

  return poly;
}

void GB2KineticIntegral::add(double coeff, double alpha0, double alpha1, const double *scales0,
                             const double *scales1) {
  double pre, gamma_inv, poly, fx0, fy0, fz0;
  double gpt_center[3], pa[3], pb[3];

  gamma_inv = 1.0 / (alpha0 + alpha1);
  pre = coeff * exp(-alpha0 * alpha1 * gamma_inv * dist_sq(r0, r1));
  compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
  pa[0] = gpt_center[0] - r0[0];
  pa[1] = gpt_center[1] - r0[1];
  pa[2] = gpt_center[2] - r0[2];
  pb[0] = gpt_center[0] - r1[0];
  pb[1] = gpt_center[1] - r1[1];
  pb[2] = gpt_center[2] - r1[2];

  i2p.reset(abs(shell_type0), abs(shell_type1));
  do {
    fx0 = gb_overlap_int1d(i2p.n0[0], i2p.n1[0], pa[0], pb[0], gamma_inv);
    fy0 = gb_overlap_int1d(i2p.n0[1], i2p.n1[1], pa[1], pb[1], gamma_inv);
    fz0 = gb_overlap_int1d(i2p.n0[2], i2p.n1[2], pa[2], pb[2], gamma_inv);
    poly = fy0 * fz0 *
        kinetic_helper(alpha0, alpha1, i2p.n0[0], i2p.n1[0], pa[0], pb[0], gamma_inv) +
        fz0 * fx0 *
            kinetic_helper(alpha0, alpha1, i2p.n0[1], i2p.n1[1], pa[1], pb[1], gamma_inv) +
        fx0 * fy0 *
            kinetic_helper(alpha0, alpha1, i2p.n0[2], i2p.n1[2], pa[2], pb[2], gamma_inv);
    work_cart[i2p.offset] += pre * scales0[i2p.ibasis0] * scales1[i2p.ibasis1] * poly;
  } while (i2p.inc());
}

/*

   GB2NuclearAttractionIntegrals

*/


GB2AttractionIntegral::GB2AttractionIntegral(long max_shell_type, double *charges,
                                             double *centers, long ncharge) :
    GB2Integral(max_shell_type), charges(charges), centers(centers), ncharge(ncharge) {
  work_g0 = new double[2 * max_shell_type + 1];
  work_g1 = new double[2 * max_shell_type + 1];
  work_g2 = new double[2 * max_shell_type + 1];
  work_boys = new double[2 * max_shell_type + 1];
}

GB2AttractionIntegral::~GB2AttractionIntegral() {
  delete[] work_g0;
  delete[] work_g1;
  delete[] work_g2;
  delete[] work_boys;
}

void
GB2AttractionIntegral::add(double coeff, double alpha0, double alpha1, const double *scales0,
                           const double *scales1) {
  double pre, gamma, gamma_inv, arg;
  double gpt_center[3], pa[3], pb[3], pc[3];

  gamma = alpha0 + alpha1;
  gamma_inv = 1.0 / gamma;
  pre = 2 * M_PI * gamma_inv * coeff * exp(-alpha0 * alpha1 * gamma_inv * dist_sq(r0, r1));
  compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
  pa[0] = gpt_center[0] - r0[0];
  pa[1] = gpt_center[1] - r0[1];
  pa[2] = gpt_center[2] - r0[2];
  pb[0] = gpt_center[0] - r1[0];
  pb[1] = gpt_center[1] - r1[1];
  pb[2] = gpt_center[2] - r1[2];

  for (long icharge = 0; icharge < ncharge; icharge++) {
    // three center for the current charge
    pc[0] = gpt_center[0] - centers[icharge * 3];
    pc[1] = gpt_center[1] - centers[icharge * 3 + 1];
    pc[2] = gpt_center[2] - centers[icharge * 3 + 2];

    // Fill the work array with the Boys function values
    arg = gamma * (pc[0] * pc[0] + pc[1] * pc[1] + pc[2] * pc[2]);
    // for (long nu=abs(shell_type0)+abs(shell_type1); nu>=0; nu--) {
    //    work_boys[nu] = boys_function(nu, arg);
    // }
    /*
       Laplace transform of the potential
    */

    int mmax = abs(shell_type0) + abs(shell_type1);
    laplace_of_potential(gamma, arg, mmax, work_boys);

    // Iterate over all combinations of Cartesian exponents
    i2p.reset(abs(shell_type0), abs(shell_type1));
    do {
      // Fill the work arrays with the polynomials
      nuclear_attraction_helper(work_g0, i2p.n0[0], i2p.n1[0], pa[0], pb[0], pc[0],
                                gamma_inv);
      nuclear_attraction_helper(work_g1, i2p.n0[1], i2p.n1[1], pa[1], pb[1], pc[1],
                                gamma_inv);
      nuclear_attraction_helper(work_g2, i2p.n0[2], i2p.n1[2], pa[2], pb[2], pc[2],
                                gamma_inv);

      // Take the product
      arg = 0;
      for (long i0 = i2p.n0[0] + i2p.n1[0]; i0 >= 0; i0--)
        for (long i1 = i2p.n0[1] + i2p.n1[1]; i1 >= 0; i1--)
          for (long i2 = i2p.n0[2] + i2p.n1[2]; i2 >= 0; i2--)
            arg += work_g0[i0] * work_g1[i1] * work_g2[i2] * work_boys[i0 + i1 + i2];

      // Finally add to the work array, accounting for opposite charge of electron and nucleus
      work_cart[i2p.offset] -=
          pre * scales0[i2p.ibasis0] * scales1[i2p.ibasis1] * arg * charges[icharge];
    } while (i2p.inc());
  }
}

void GB2NuclearAttractionIntegral::laplace_of_potential(double gamma, double arg, long mmax,
                                                        double *output) {
  boys_function_array(mmax, arg, output);
}

void GB2ErfAttractionIntegral::laplace_of_potential(double gamma, double arg, long mmax,
                                                    double *output) {
  double efac = mu * mu / (mu * mu + gamma);
  boys_function_array(mmax, arg * efac, output);
  double prefac = sqrt(efac);
  for (long m = 0; m <= mmax; m++) {
    output[m] *= prefac;
    prefac *= efac;
  }
}

void GB2GaussAttractionIntegral::laplace_of_potential(double gamma, double arg, long mmax,
                                                      double *output) {
  double afac = alpha / (gamma + alpha);
  double prefac =
      (M_PI / (gamma + alpha)) * sqrt(M_PI / (gamma + alpha)) * c * exp(-arg * afac) *
          (gamma / (2 * M_PI));
  for (long m = 0; m <= mmax; m++) {
    output[m] = prefac;
    prefac *= afac;
  }
}

/*

 GB2MomentIntegral

 */


double moment_helper(long n0, long n1, long n2, double pa, double pb, double pc, double gamma_inv) {
  /*
   The Obara Saika Scheme, equations A7 and A8 in "Efficient recursive computation of
   molecular integrals over Cartesian Gaussian functions", S. Obara, A. Saika, Journal
   of Chemical Physics, vol. 84, p. 3963, y. 1986.

   Note that we use a shift in the indices:
   (0_A|R(0)|0_B) in the paper is stored in work_mm[1][1][1].
   In this allows to set all terms in with the angular momentum index becomes -1 to 0,
   avoiding a lot of if-statements.
  */

  double result;
  // if n2 == 0, we just need the overlap.
  if (n2 == 0) {
    result = gb_overlap_int1d(n0, n1, pa, pb, (2.0e0 * gamma_inv));
  } else {
    long m, l, n;
    double work_mm[n0 + 2][n1 + 2][n2 + 2];

    for (l = 0; l < (n0 + 2); l++) {
      for (m = 0; m < (n1 + 2); m++) {
        for (n = 0; n < (n2 + 2); n++) {
          work_mm[l][m][n] = 0.0e0;
        }
      }
    }

    // The auxiliary overlap (0_A|R(0)|0_B) is stored in work_mm[1][1][1]:
    work_mm[1][1][1] = gb_overlap_int1d(0, 0, pa, pb, (2.0e0 * gamma_inv));

    // Equation A8 in the Obara-Saika paper for (0_A|R(mu + 1)|0_B):
    for (n = 1; n < (n2 + 1); n++) {
      work_mm[1][1][n + 1] =
          pc * work_mm[1][1][n] + (gamma_inv * (n - 1) * work_mm[1][1][n - 1]);
    }

    // Equation A7 in the Obara-Saika paper for (0_A|R(mu)|b + 1):
    for (m = 1; m < (n1 + 1); m++) {
      for (n = 1; n <= (n2 + 1); n++) {
        work_mm[1][m + 1][n] = pb * work_mm[1][m][n]
            + (gamma_inv * (m - 1) * work_mm[1][m - 1][n])
            + (gamma_inv * (n - 1) * work_mm[1][m][n - 1]);
      }
    }

    // Equation A7 in the Obara-Saika paper for (a + 1|R(mu)|b):
    for (l = 1; l < (n0 + 1); l++) {
      for (m = 1; m <= (n1 + 1); m++) {
        for (n = 1; n <= (n2 + 1); n++) {
          work_mm[l + 1][m][n] = pa * work_mm[l][m][n]
              + (gamma_inv * (l - 1) * work_mm[l - 1][m][n])
              + (gamma_inv * (m - 1) * work_mm[l][m - 1][n])
              + (gamma_inv * (n - 1) * work_mm[l][m][n - 1]);
        }
      }
    }

    result = work_mm[n0 + 1][n1 + 1][n2 + 1];
  }

  return result;
}

GB2MomentIntegral::GB2MomentIntegral(long max_shell_type, long *xyz, double *center)
    : GB2Integral(max_shell_type), xyz(xyz), center(center) {
  if (xyz[0] + xyz[1] + xyz[2] < 0)
    throw std::domain_error(" sum < 0");
  if ((xyz[0] < 0) || (xyz[1] < 0) || (xyz[2] < 0))
    throw std::domain_error(" all elements of xyz must be >= 0");
}

void GB2MomentIntegral::add(double coeff, double alpha0, double alpha1,
                            const double *scales0, const double *scales1) {
  double pre, gamma_inv, twogamma_inv;
  double gpt_center[3], pa[3], pb[3], pc[3];

  gamma_inv = 1.0 / (alpha0 + alpha1);
  twogamma_inv = 0.50 / (alpha0 + alpha1);
  pre = coeff * exp(-alpha0 * alpha1 * gamma_inv * dist_sq(r0, r1));
  compute_gpt_center(alpha0, r0, alpha1, r1, gamma_inv, gpt_center);
  pa[0] = gpt_center[0] - r0[0];
  pa[1] = gpt_center[1] - r0[1];
  pa[2] = gpt_center[2] - r0[2];
  pb[0] = gpt_center[0] - r1[0];
  pb[1] = gpt_center[1] - r1[1];
  pb[2] = gpt_center[2] - r1[2];
  pc[0] = gpt_center[0] - center[0];
  pc[1] = gpt_center[1] - center[1];
  pc[2] = gpt_center[2] - center[2];
  i2p.reset(abs(shell_type0), abs(shell_type1));

  do {
    work_cart[i2p.offset] += pre * (
        moment_helper(i2p.n0[0], i2p.n1[0], xyz[0], pa[0], pb[0], pc[0], twogamma_inv) *
            moment_helper(i2p.n0[1], i2p.n1[1], xyz[1], pa[1], pb[1], pc[1], twogamma_inv) *
            moment_helper(i2p.n0[2], i2p.n1[2], xyz[2], pa[2], pb[2], pc[2], twogamma_inv) *
            scales0[i2p.ibasis0] * scales1[i2p.ibasis1]);
  } while (i2p.inc());
}
