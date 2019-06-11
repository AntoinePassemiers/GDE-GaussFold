# -*- coding: utf-8 -*-
# chain.py: Base class for amino acid chains
# author : Antoine Passemiers

from gaussfold.aa import *


class Chain:

    def __init__(self):
        self.amino_acids = list()

    def add(self, amino_acid):
        self.amino_acids.append(amino_acid)
        # TODO: add peptide bonds

    @staticmethod
    def from_string(s):
        abb_to_clk = dict()
        for ckl in AminoAcid.__subclasses__():
            abb_to_clk[ckl().abbreviation] = ckl

        chain = Chain()
        for c in s:
            chain.add(abb_to_clk[c])
        return chain
