# -*- coding: utf-8 -*-
# amino_acid.py: Base class for amino acids
# author : Antoine Passemiers

from gaussfold.atom import Atom, Bond
from gaussfold.atom import Carbon, Hydrogen, Oxygen, Nitrogen

from abc import ABCMeta, abstractmethod


class AminoAcid:

    def __init__(self, aa_name, abbreviation, explicit_hydrogens=False):
        self.aa_name = aa_name
        self.abbreviation = abbreviation
        self.explicit_hydrogens = explicit_hydrogens

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

    def add_atom(self, atom):
        assert(isinstance(atom, Atom))
        if atom not in self.atoms:
            if (not isinstance(atom, Hydrogen)) or self.explicit_hydrogens:
                self.atoms.append(atom)

    def add_bond(self, bond):
        assert(isinstance(bond, Bond))
        is_1_h = isinstance(bond.atom1, Hydrogen)
        is_2_h = isinstance(bond.atom2, Hydrogen)
        if (not is_1_h) or (not is_2_h) or self.explicit_hydrogens:
            self.bonds.append(bond)

    def __to_pdb__(self, serial, chain_id, res_seq):
        s = ''
        for atom in self.atoms:
            s += atom.__to_pdb__(serial, chain_id, res_seq)
            serial += 1
        return s, serial
