# -*- coding: utf-8 -*-
# chain.py: Base class for amino acid chains
# author : Antoine Passemiers

from gaussfold.aa import *
from gaussfold.atom import Bond


class Chain:

    def __init__(self):
        self.amino_acids = list()

    def add(self, amino_acid):
        if len(self.amino_acids) > 0:
            peptide_bond = Bond(self.amino_acids[-1].C, amino_acid.N)
            self.amino_acids[-1].add_bond(peptide_bond)
        self.amino_acids.append(amino_acid)

    def atoms(self):
        atom_list = list()
        for amino_acid in self.amino_acids:
            atom_list += list(amino_acid.atoms)
        return atom_list

    @staticmethod
    def from_string(s):
        abb_to_clk = dict()
        for ckl in AminoAcid.__subclasses__():
            abb_to_clk[ckl().abbreviation] = ckl

        chain = Chain()
        for c in s:
            chain.add(abb_to_clk[c]())
        return chain
