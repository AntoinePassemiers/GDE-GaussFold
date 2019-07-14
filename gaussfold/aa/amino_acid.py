# -*- coding: utf-8 -*-
# amino_acid.py: Base class for amino acids
# author : Antoine Passemiers

from gaussfold.atom import Atom, Bond
from gaussfold.atom import Carbon, Hydrogen, Oxygen, Nitrogen

import numpy as np
from abc import ABCMeta, abstractmethod


class AminoAcid:

    def __init__(self, aa_name, abbreviation, explicit_hydrogens=False, c='CA'):
        self.aa_name = aa_name
        self.abbreviation = abbreviation
        self.explicit_hydrogens = explicit_hydrogens
        self.c = c
        assert(self.c in ['CA', 'CB'])

        self.atoms = list()
        self.bonds = list()

        # Define backbone atoms

        self.N = Nitrogen('N')
        self.add_atom(self.N)

        self.CA = Carbon('CA')
        self.add_atom(self.CA)

        self.C = Carbon('C')
        self.add_atom(self.C)

        self.O = Oxygen('O')
        self.add_atom(self.O)

        self.add_bond(Bond(self.N, self.CA))
        self.add_bond(Bond(self.CA, self.C))
        self.add_bond(Bond(self.C, self.O, order=2))

        # Add hydrogen

        self.H = Hydrogen('H')
        self.add_atom(self.H)
        self.add_bond(Bond(self.H, self.CA))

    def ref(self):
        if self.c == 'CB' and hasattr(self, 'CB'):
            return self.CB
        else:
            return self.CA

    def phi(self, left_atom):
        return self.dihedral_angle(
                left_atom.get_coords(),
                self.N.get_coords(),
                self.CA.get_coords(),
                self.C.get_coords())

    def psi(self, right_atom):
        return self.dihedral_angle(
                self.N.get_coords(),
                self.CA.get_coords(),
                self.C.get_coords(),
                right_atom.get_coords())

    def dihedral_angle(self, p1, p2, p3, p4):
        
        # Compute vectors from coordinates
        v1 = p2 - p1
        v1 /= np.linalg.norm(v1)
        v2 = p3 - p2
        v2 /= np.linalg.norm(v2)
        v3 = p4 - p3
        v3 /= np.linalg.norm(v3)

        # Normal vector of the plane containing v1 and v2
        nv1 = np.cross(v1, v2)

        # Normal vector of the plane containing v2 and v3
        nv2 = np.cross(v2, v3)

        # Vector that is orthogonal to nv1 and v2
        o1 = np.cross(nv1, v2)

        # Coordinates of nv2 in the frame (nv1, v2, o1)
        x = np.dot(nv1, nv2)
        y = np.dot(o1, nv2)

        return np.arctan2(y, x) * 180. / np.pi

    def add_atom(self, atom):
        assert(isinstance(atom, Atom))
        if atom not in self.atoms:
            if (not isinstance(atom, Hydrogen)) or self.explicit_hydrogens:
                self.atoms.append(atom)

    def add_bond(self, bond):
        assert(isinstance(bond, Bond))
        is_1_h = isinstance(bond.atom1, Hydrogen)
        is_2_h = isinstance(bond.atom2, Hydrogen)
        if ((not is_1_h) and (not is_2_h)) or self.explicit_hydrogens:
            self.bonds.append(bond)

    def __to_pdb__(self, serial, chain_id, res_seq):
        s = ''
        for atom in self.atoms:
            s += atom.__to_pdb__(serial, chain_id, res_seq)
            serial += 1
        return s, serial
