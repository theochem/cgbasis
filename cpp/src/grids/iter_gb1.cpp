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


#include <cstdlib>
#include <cstring>
#include "../common.h"
#include "iter_gb1.h"

IterGB1::IterGB1(GBasis *gbasis)
    : gbasis(gbasis), basis_offsets(gbasis->get_basis_offsets()),
      // public fields
      shell_type0(0), con_coeff(0.0), alpha0(0.0), r0(NULL), scales0(NULL),
      ibasis0(0),
      // internal fields
      ishell0(0), nprim0(0), iprim0(0), oprim0(0) {}

int IterGB1::inc_shell() {
  // Increment shell and related counters.
  if (ishell0 < gbasis->nshell - 1) {
    oprim0 += nprim0;
    ishell0++;
    update_shell();
    return 1;
  } else {
    ishell0 = 0;
    oprim0 = 0;
    update_shell();
    return 0;
  }
}

void IterGB1::update_shell() {
  // Update fields that depend on shell and related counters.
  nprim0 = gbasis->nprims[ishell0];
  // Update indexes in output array
  ibasis0 = basis_offsets[ishell0];
  // update centers
  r0 = gbasis->centers + 3 * gbasis->shell_map[ishell0];
  // update shell types
  shell_type0 = gbasis->shell_types[ishell0];
  // reset contraction counters
  iprim0 = 0;
}

int IterGB1::inc_prim() {
  // Increment primitive counters.
  if (iprim0 < nprim0 - 1) {
    iprim0++;
    update_prim();
    return 1;
  } else {
    iprim0 = 0;
    update_prim();
    return 0;
  }
}

void IterGB1::update_prim() {
  // Update fields that depend on primitive counters.
  alpha0 = gbasis->alphas[oprim0 + iprim0];
  con_coeff = gbasis->con_coeffs[oprim0 + iprim0];
  scales0 = gbasis->get_scales(oprim0 + iprim0);
}

void IterGB1::store(const double *work, double *output, long dim) {
  // This routine is hardwired to work only for the dense storage
  const long n0 = get_shell_nbasis(shell_type0);
  memcpy(output + ibasis0 * dim, work, n0 * dim * sizeof(double));
}
