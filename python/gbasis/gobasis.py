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
"""Gaussian orbital basis set module."""
from __future__ import annotations

from typing import Union, Dict, Iterable, List

import numpy as np

from .cext import GOBasis
from .iobas import load_basis_atom_map_nwchem, load_basis_atom_map_gbs, dump_basis_atom_map_gbs
from .periodic import sym2num
from .utils import typecheck_geo, to_bset_path

__all__ = [
    'get_gobasis', 'GOBasisDesc', 'GOBasisFamily', 'go_basis_families',
    'GOBasisAtom', 'GOBasisContraction',
]


def get_gobasis(coordinates: np.ndarray, numbers: np.ndarray, default: Union[str, GOBasisFamily],
                element_map: Dict = None, index_map: Dict = None, pure: bool = True) -> GOBasis:
    """Return GOBasis for a given molecule. This is the standard way to get integrals.

    Parameters
    ----------
    coordinates
        A (N, 3) float numpy array with Cartesian coordinates of the atoms in Bohr.
    numbers
        A (N,) numpy vector with the atomic numbers.
    default
        The default basis set applied to each atom.
    element_map
        A dictionary with element names or numbers as keys, and basis sets
        as values. These specs override the default basis.
    index_map
        A dictionary with atomic indexes (based on the order of the atoms)
        as keys and basis sets as values.
    pure
        By default pure basis functions are used. Set this to false to
        switch to Cartesian basis functions.

    Returns
    -------
        A GOBasis instance which can be used to calculate integrals for the molecule.

    Examples
    --------
    >>>> gbasis = get_gobasis(np.array([[0,0,0], [0,0,1]]), np.array([1,1]), "sto-3g")

    >>>> gbasis.compute_overlap()

    """
    gobasis_desc = GOBasisDesc(default, element_map, index_map, pure)
    return gobasis_desc.apply_to(coordinates, numbers)


class GOBasisDesc:
    """A user specification of the basis set."""

    def __init__(self, default: Union[str, GOBasisFamily],
                 element_map: Dict = None, index_map: Dict = None,
                 pure: bool = True):
        """

        Parameters
        ----------
        default
            The default basis set applied to each atom.
        element_map
            A dictionary with element names or numbers as keys, and
            basis sets as values. These specs override the default basis
        index_map
            A dictionary with atomic indexes (based on the order of the
            atoms) as keys and basis sets as values
        pure
            By default pure basis functions are used. Set this to false to
            switch to Cartesian basis functions.
        """
        self.default = default
        if element_map is None:
            self.element_map = {}
        else:
            self.element_map = element_map
        if index_map is None:
            self.index_map = {}
        else:
            self.index_map = index_map
        self.pure = pure

        # Update the element map such that it only contains numbers as keys.
        for key in list(self.element_map.keys()):
            if isinstance(key, str):
                number = sym2num[key]
                self.element_map[number] = element_map[key]
                del element_map[key]

    def apply_to(self, coordinates: np.ndarray, numbers: np.ndarray) -> GOBasis:
        """Construct a GOBasis object for the given molecular geometry

        Parameters
        ----------
        coordinates
            A (N, 3) float numpy array with Cartesian coordinates of the
            atoms in Bohr.
        numbers
            A (N,) numpy vector with the atomic numbers.

        Returns
        -------
        A GOBasis object.

        Note that the geometry specified by the arguments may also contain
        ghost atoms.

        """
        natom, coordinates, numbers = typecheck_geo(coordinates, numbers, need_pseudo_numbers=False)

        shell_map = []
        nprims = []
        shell_types = []
        alphas = []
        con_coeffs = []

        def get_basis(i: int, n: int) -> str:
            """Look up the basis for a given atom. First by `index_map`, then by `element_map`,
            then by `default`.

            Parameters
            ----------
            i
                The index of the atom.
            n
                The atom number of the atom.

            Returns
            -------
            The basis set of that atom.

            Raises
            ------
            KeyError
                If the basis set doesn't contain the atom, and default isn't provided.

            """
            basis = self.index_map.get(i)
            if basis is not None:
                return basis
            basis = self.element_map.get(n)
            if basis is not None:
                return basis
            if self.default is None:
                raise KeyError(f'Could not find basis for atom {i}.')
            else:
                return self.default

        def translate_basis(basis_x: Union[str, GOBasisFamily, GOBasisAtom], n: int) -> GOBasisAtom:
            """
            Get a GOBasisAtom instance of an atomic element.

            Parameters
            ----------
            basis_x
                The basis specification.
            n
                The atomic number of the atom in question

            Returns
            -------
            An instance of GOBasisAtom.

            Raises
            ------
            ValueError
                If the basis set doesn't exist or doesn't contain the atom.

            """
            if isinstance(basis_x, str):
                basis_fam = go_basis_families.get(basis_x.lower())
                if basis_fam is None:
                    raise ValueError(f'Unknown basis family: {basis_x}')
                basis_atom = basis_fam.get(n)
                if basis_atom is None:
                    raise ValueError(
                        f'The basis family {basis_x} does not contain element {n}.')
                return basis_atom
            elif isinstance(basis_x, GOBasisFamily):
                basis_atom = basis_x.get(n)
                if basis_atom is None:
                    raise ValueError(
                        f'The basis family {basis_x.name} does not contain element {n}.')
                return basis_atom
            elif isinstance(basis_x, GOBasisAtom):
                return basis_x
            else:
                raise ValueError(f'Can not interpret {basis_x} as an atomic basis function.')

        # Loop over the atoms and fill in all the lists
        for i in range(natom):
            n = numbers[i]
            basis_x = get_basis(i, n)
            basis_atom = translate_basis(basis_x, n)
            basis_atom.extend(i, shell_map, nprims, shell_types, alphas, con_coeffs, self.pure)

        # Return the Gaussian basis object.
        return GOBasis(coordinates, shell_map, nprims, shell_types, alphas, con_coeffs)


class GOBasisFamily:
    """An object to store a basis set.

    Attributes
    ----------
    name
        The name of the basis set. It will be printed in output later.
    basis_atom_map
        The dictionary of atomic numbers as keys and GOBasisAtom instances as values.
    filename
        The filename storing the basis set.
    """

    def __init__(self, name: str, basis_atom_map: Dict[int, GOBasisAtom] = None,
                 filename: str = None):
        """
        Parameters
        ----------
        name
            The name of the basis set family, e.g. 'STO-3G'.
        basis_atom_map
            A dictionary with atomic numbers as keys and GOBasisAtom
            instances as values. This can be used to craft your own basis set.
        filename
            A file to load the basis set from when needed. Useful for specifying
            a basis set that doesn't already exist in the library. The extension of
            the file must be `.nwchem` or `.gbs`.

        Raises
        ------
        ValueError
            If neither `basis_atom_map` nor `filename` are provided.

        Examples
        --------
        >>>> mybasisset = GOBasisFamily("complete_basis", filename="complete_basis.nwchem")

        >>>> get_gobasis(np.array([[0,0,0],]), np.array([1]), default=mybasisset)

        """
        if basis_atom_map is None and filename is None:
            raise ValueError('Either one of the two optional arguments must be provided.')
        self.name = name
        self.basis_atom_map = basis_atom_map
        self.filename = filename

    def get(self, number: int) -> GOBasisAtom:
        """Get the GOBasisAtom instance for an atomic number.

        Parameters
        ----------
        number
            The atomic number

        Returns
        -------
            An instance of GOBasisAtom
        """
        if self.basis_atom_map is None:
            self.load()
        return self.basis_atom_map[number]

    def load(self):
        """Load the basis set from file if it hasn't been done already.
        If the basis_atom_map is already defined (not None), then the load method is ignored.
        """
        # if basis_atom_map is already defined
        if self.basis_atom_map is not None:
            # ignore load method
            return
        if self.filename.endswith('.nwchem'):
            self.basis_atom_map = load_basis_atom_map_nwchem(self.filename)
        elif self.filename.endswith('.gbs'):
            self.basis_atom_map = load_basis_atom_map_gbs(self.filename)
        else:
            raise IOError(f'File format not supported: {self.filename}')
        self._to_arrays()
        self._to_segmented()
        self._normalize_contractions()

    def dump(self, filename: str):
        """Dump the basis set in the gbs format.

        Parameters
        ----------
        filename
            Name of the gbs file that will be created.
        """
        self.load()
        if filename.endswith('.gbs'):
            dump_basis_atom_map_gbs(filename, self.name, self.basis_atom_map)
        else:
            raise IOError(f'File format not supported: {filename}')

    def _to_arrays(self):
        """Convert all contraction attributes to numpy arrays."""
        for ba in self.basis_atom_map.values():
            for bc in ba.bcs:
                bc.to_arrays()

    def _to_segmented(self):
        """Convert all contractions from generalized to segmented"""
        new_basis_atom_map = {}
        for n, ba in self.basis_atom_map.items():
            new_bcs = []
            for bc in ba.bcs:
                new_bcs.extend(bc.get_segmented_bcs())
            new_ba = GOBasisAtom(new_bcs)
            new_basis_atom_map[n] = new_ba
        self.basis_atom_map = new_basis_atom_map

    def _normalize_contractions(self):
        """Renormalize all contractions."""
        for ba in self.basis_atom_map.values():
            for bc in ba.bcs:
                bc.normalize()


go_basis_families_list = [
    GOBasisFamily('STO-3G', filename=to_bset_path('sto-3g.nwchem')),
    GOBasisFamily('STO-6G', filename=to_bset_path('sto-6g.nwchem')),
    GOBasisFamily('3-21G', filename=to_bset_path('3-21g.nwchem')),
    GOBasisFamily('3-21G(d)', filename=to_bset_path('3-21g(d).nwchem')),
    GOBasisFamily('3-21++G(d)', filename=to_bset_path('3-21++g(d).nwchem')),
    GOBasisFamily('4-31G', filename=to_bset_path('4-31g.nwchem')),
    GOBasisFamily('6-31G', filename=to_bset_path('6-31g.nwchem')),
    GOBasisFamily('6-31G(d)', filename=to_bset_path('6-31g(d).nwchem')),
    GOBasisFamily('6-31G(d,p)', filename=to_bset_path('6-31g(d,p).nwchem')),
    GOBasisFamily('6-31+G', filename=to_bset_path('6-31+g.nwchem')),
    GOBasisFamily('6-31+G(d)', filename=to_bset_path('6-31+g(d).nwchem')),
    GOBasisFamily('6-31++G(d,p)', filename=to_bset_path('6-31++g(d,p).nwchem')),
    GOBasisFamily('cc-pVDZ', filename=to_bset_path('cc-pvdz.nwchem')),
    GOBasisFamily('cc-pVTZ', filename=to_bset_path('cc-pvtz.nwchem')),
    GOBasisFamily('cc-pVQZ', filename=to_bset_path('cc-pvqz.nwchem')),
    GOBasisFamily('cc-pCVDZ', filename=to_bset_path('cc-pcvdz.nwchem')),
    GOBasisFamily('cc-pCVTZ', filename=to_bset_path('cc-pcvtz.nwchem')),
    GOBasisFamily('cc-pCVQZ', filename=to_bset_path('cc-pcvqz.nwchem')),
    GOBasisFamily('aug-cc-pVDZ', filename=to_bset_path('aug-cc-pvdz.nwchem')),
    GOBasisFamily('aug-cc-pVTZ', filename=to_bset_path('aug-cc-pvtz.nwchem')),
    GOBasisFamily('aug-cc-pVQZ', filename=to_bset_path('aug-cc-pvqz.nwchem')),
    GOBasisFamily('aug-cc-pV5Z', filename=to_bset_path('aug-cc-pv5z.nwchem')),
    GOBasisFamily('aug-cc-pV6Z', filename=to_bset_path('aug-cc-pv6z.nwchem')),
    GOBasisFamily('aug-cc-pCVDZ', filename=to_bset_path('aug-cc-pcvdz.nwchem')),
    GOBasisFamily('aug-cc-pCVTZ', filename=to_bset_path('aug-cc-pcvtz.nwchem')),
    GOBasisFamily('aug-cc-pCVQZ', filename=to_bset_path('aug-cc-pcvqz.nwchem')),
    GOBasisFamily('def2-svpd', filename=to_bset_path('def2-svpd.nwchem')),
    GOBasisFamily('def2-tzvp', filename=to_bset_path('def2-tzvp.nwchem')),
    GOBasisFamily('def2-tzvpd', filename=to_bset_path('def2-tzvpd.nwchem')),
    GOBasisFamily('def2-qzvp', filename=to_bset_path('def2-qzvp.nwchem')),
    GOBasisFamily('def2-qzvpd', filename=to_bset_path('def2-qzvpd.nwchem')),
    GOBasisFamily('ANO-RCC', filename=to_bset_path('ano-rcc.nwchem')),
]
go_basis_families = dict((bf.name.lower(), bf) for bf in go_basis_families_list)


class GOBasisAtom:
    """Description of an atomic basis set with segmented contractions."""

    def __init__(self, bcs: Iterable):
        """
        Parameters
        ----------
        bcs
            GOBasisContraction instances

        """
        self.bcs = bcs

    def extend(self, i: int, shell_map: List[int], nprims: List[int], shell_types: List[int],
               alphas: List[float], con_coeffs: List[float], pure=True):
        """
        Add basis functions to an atom. This can take an existing set of parameters for GOBasis and
        add an atom (which is defined by this instance) to it.

        Parameters
        ----------
        i
            The index of the center of this atom.
        shell_map
            The index of the centres of each shell. This instance's centre will be appended to it.
        nprims
            The number of primitives in each shell. This instance's centre will be appended to it.
        shell_types
            The angular momentum of the shell. 0 = S, 1 = P, 2 = D (cartesian), -2 = D (pure)
            3 = F (cartesian), -3 = F (pure). This instance's shells will be appended to it.
        alphas
            The exponent of each primitive.
        con_coeffs
            The contraction coefficient of each primitive.
        pure
            By default pure basis functions are used. Set this to False to
            switch to Cartesian basis functions.

        """
        for bc in self.bcs:
            if bc.is_generalized():
                raise ValueError('Generalized contractions are not supported (yet).')
            shell_map.append(i)
            nprims.append(len(bc.alphas))
            if pure and bc.shell_type >= 2:
                shell_types.append(-bc.shell_type)
            else:
                shell_types.append(bc.shell_type)
            alphas.extend(bc.alphas)
            con_coeffs.extend(bc.con_coeffs)


class GOBasisContraction:
    """A basis set shell. Can represent either a segmented or a generalized contraction."""

    def __init__(self, shell_type: int, alphas: Iterable, con_coeffs: Iterable):
        """
        Parameters
        ----------
        shell_type
            The angular momentum quantum number of the shell. (0=s, 1=p, 2=d(cartesian),
            -2=d(pure)).
        alphas
            The exponents of each primitive. Should be 1D array-like object.
        con_coeffs
            The contraction coefficients for the shell.

            In the case of a segmented basis set, this is just a 1D array-like object
            with the same size as alphas.

            In the case of a generalized contraction, this is a
            2D array-like object, where each row corresponds to a primitive and
            the columns correspond to different contractions.
        """
        self.shell_type = shell_type
        self.alphas = alphas
        self.con_coeffs = con_coeffs

    def to_arrays(self):
        """Convert the alphas and con_coeffs attributes to numpy arrays."""
        self.alphas = np.array(self.alphas)
        self.con_coeffs = np.array(self.con_coeffs)

    def is_generalized(self) -> bool:
        """Return True if this is a generalized contraction."""
        return self.con_coeffs.ndim >= 2

    def get_segmented_bcs(self) -> List[GOBasisContraction]:
        """Return a list of segmented contractions if the original instance contains a generalized
        contraction.

        Raises
        ------
        TypeError
            If contraction is not generalized.
        """
        if not self.is_generalized():
            raise TypeError('Conversion to segmented contractions only makes sense for '
                            'generalized contractions.')
        return [
            GOBasisContraction(self.shell_type, self.alphas, self.con_coeffs[:, i])
            for i in range(self.con_coeffs.shape[1])
        ]

    def normalize(self):
        """Normalize the contraction.

        Raises
        ------
        NotImplementedError
            If the contraction is not segmented.
        """
        if self.is_generalized():
            raise NotImplementedError("Only segmented contractions can be normalized.")
        # Warning! Ugly code ahead to avoid re-implementing the norm of contraction. The
        # code below (ab)uses the GOBasis machinery to get that result.
        # 1) Construct a GOBasis object with only this contraction.
        centers = np.array([[0.0, 0.0, 0.0]])
        shell_map = np.array([0])
        nprims = np.array([len(self.alphas)])
        shell_types = np.array([self.shell_type])
        alphas = self.alphas
        con_coeffs = self.con_coeffs
        gobasis = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        # 2) Get the first diagonal element of the overlap matrix
        olpdiag = gobasis.compute_overlap()[0, 0]
        # 3) Normalize the contraction
        self.con_coeffs /= np.sqrt(olpdiag)
