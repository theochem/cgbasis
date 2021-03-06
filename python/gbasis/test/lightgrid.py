# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
"""Simple molecule integration grids."""

from __future__ import annotations

from functools import reduce
from typing import Union, Tuple

import numpy as np

__all__ = ['generate_molecular_grid', 'integrate']


class Shell:
    """A 1s-type Slater density shell."""

    def __init__(self, population: float, exponent: float, center: np.ndarray = None):
        """Initialize a Shell instance.

        Parameters
        ----------
        population
            The number of electrons in the shell.
        exponent
            The exponent of the Slater function.
        center
            The atomic center for the shell. An array of (3, ) shape.
        """
        self.population = population
        self.exponent = exponent
        self.center = center

    def clone(self, center: np.ndarray) -> Shell:
        """Make a copy of the shell with a new center.
        Parameters
        ----------
        center
            The atomic center for the shell. An array of (3, ) shape.

        Returns
        -------
        Another instance of Shell with the new center.
        """
        return Shell(self.population, self.exponent, center)

    def compute_density(self, distance: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Compute the density of the shell as function of the distance from the centre.

        Parameters
        ----------
        distance
            The distance from the centre.

        Returns
        -------

        """
        prefactor = self.population * self.exponent ** 3 / (8 * np.pi)
        return prefactor * np.exp(-distance * self.exponent)

    def get_random_points(self, npoint: int) -> np.ndarray:
        """Sample npoint random points from shape function of the shell."""
        radii = np.random.gamma(3, 1.0 / self.exponent, npoint)
        points = np.zeros((npoint, 3))
        norms = np.zeros(npoint)
        while True:
            mask = norms < 0.1
            nnew = mask.sum()
            if nnew == 0:
                break
            newpoints = np.random.normal(0, 1, (nnew, 3))
            newnorms = np.sqrt(np.einsum('ij,ij->i', newpoints, newpoints))
            points[mask] = newpoints
            norms[mask] = newnorms
        return points * (radii / norms).reshape(-1, 1) + self.center


setups = {
    1: [Shell(1.0, 2.0)],
    3: [Shell(1.86359, 5.56763), Shell(1.13641, 0.80520)],
    6: [Shell(1.70730, 12.79758), Shell(4.29270, 1.85580)],
    7: [Shell(1.68283, 15.13096), Shell(5.31717, 2.19942)],
    8: [Shell(1.66122, 17.46129), Shell(6.33878, 2.54326)],
    9: [Shell(1.64171, 19.78991), Shell(7.35829, 2.88601)],
}


def generate_molecular_grid(numbers: np.ndarray, coordinates: np.ndarray,
                            npoint_per_electron: int = 100, seed: int = 1) -> Tuple[np.ndarray,
                                                                                    np.ndarray]:
    """Generate a molecular integration grid, using importance sampling.

    Parameters
    ----------
    numbers
        Atomic numbers. shape=(natom,), dtype=int
    coordinates
        Atomic coordinates in Bohr. shape(natom, 3), dtype=float
    npoint_per_electron
        The number of grid points per atomic shell.
    seed
        The random seed to use when generating the grid.

    Returns
    -------
    points
        Positions of the grid points. shape=(npoint, 3), dtype=float
    weights
        Integration grid weights for all points. shape=(npoint,), dtype=float
    """
    state = np.random.get_state()
    np.random.seed(seed)
    try:
        # Get a list of shells
        shells = []
        for number, coordinate in zip(numbers, coordinates):
            for shell in setups[number]:
                shells.append(shell.clone(coordinate))

        # Initialize outputs
        shell_points = []
        shell_weights = []

        # Random sampling
        nelec = sum(shell.population for shell in shells)
        weights_sum = 0.0
        for ishell0, shell0 in enumerate(shells):
            npoint_shell = int(np.round(npoint_per_electron * shell0.population))
            shell_points.append(shell0.get_random_points(npoint_shell))
            prob = 0.0
            for ishell1, shell1 in enumerate(shells):
                deltas = shell_points[-1] - shell1.center
                distances = np.sqrt(np.einsum('ij,ij->i', deltas, deltas))
                prob += shell1.compute_density(distances)
            prob /= nelec
            scale = shell0.population / npoint_shell
            shell_weights.append(scale / prob)
            weights_sum += shell0.population

        # Fix shapes
        points = np.concatenate(shell_points)
        weights = np.concatenate(shell_weights)
    finally:
        np.random.set_state(state)

    return points, weights / weights_sum


def integrate(*args: np.ndarray) -> np.float64:
    """Replaces grid.integrate for gbasis tests.
    Simply takes a dot product of all 1D numpy arrays passed to it.

    Parameters
    ----------
    args
        1D numpy arrays
    Returns
    -------
        The integral of the grids.
    """
    return np.sum(reduce(np.multiply, args))
