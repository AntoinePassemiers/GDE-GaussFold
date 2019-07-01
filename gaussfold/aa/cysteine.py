# -*- coding: utf-8 -*-
# cysteine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Hydrogen, Sulfur


class Cysteine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'CYS', 'C', **kwargs)

        # Add side chain atoms

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.SG = Sulfur('SG')
        self.add_atom(self.SG)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.SG, self.CB))

        # Add hydrogens

        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CB))

        self.H2 = Hydrogen('H2')
        self.add_atom(self.H2)
        self.add_bond(Bond(self.H2, self.CB))

        self.H3 = Hydrogen('H3')
        self.add_atom(self.H3)
        self.add_bond(Bond(self.H3, self.SG))
