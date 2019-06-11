# -*- coding: utf-8 -*-
# amino_acid.py: Base class for amino acids
# author : Antoine Passemiers

from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen

from abc import ABCMeta, abstractmethod


class AminoAcid:

    def __init__(self, aa_name, abbreviation):
        self.aa_name = aa_name
        self.abbreviation = abbreviation

        self.atoms = list()
        self.bonds = list()

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

    def add_atom(self, atom):
        if atom not in self.atoms:
            self.atoms.append(atom)

    def add_bond(self, bond):
        self.bonds.append(bond)

    def __to_pdb__(self, serial, chain_id, res_seq):
        s = ''
        for atom in self.atoms:
            s += atom.__to_pdb__(serial, chain_id, res_seq)
            serial += 1
        return s, serial
