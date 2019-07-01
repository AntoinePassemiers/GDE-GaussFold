# -*- coding: utf-8 -*-
# chain.py: Base class for amino acid chains
# author : Antoine Passemiers

from gaussfold.aa import *
from gaussfold.atom import Bond


class Chain:

    def __init__(self):
        self.amino_acids = list()

    def update(self):
        identifier = 0
        for amino_acid in self.amino_acids:
            for atom in amino_acid.atoms:
                atom.identifier = identifier
                identifier += 1

    def add(self, amino_acid):
        if len(self.amino_acids) > 0:
            peptide_bond = Bond(self.amino_acids[-1].C, amino_acid.N)
            self.amino_acids[-1].add_bond(peptide_bond)
        self.amino_acids.append(amino_acid)
        self.update()

    def atoms(self):
        atom_list = list()
        for amino_acid in self.amino_acids:
            atom_list += list(amino_acid.atoms)
        return atom_list

    def bonds(self):
        bond_list = list()
        for amino_acid in self.amino_acids:
            bond_list += list(amino_acid.bonds)
        return bond_list

    @staticmethod
    def from_string(s):
        abb_to_clk = dict()
        for ckl in AminoAcid.__subclasses__():
            abb_to_clk[ckl().abbreviation] = ckl

        chain = Chain()
        for c in s:
            chain.add(abb_to_clk[c]())
        return chain

    def __getitem__(self, key):
        return self.amino_acids[key]

    def __setitem__(self, key, value):
        self.amino_acids[key] = value
        # TODO: check value type
        self.update()

    def __len__(self):
        return len(self.amino_acids)
