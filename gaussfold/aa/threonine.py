# -*- coding: utf-8 -*-
# threonine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Threonine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'THR', 'T')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.OG1 = Oxygen('OG1')
        self.add_atom(self.OG1)

        self.CG2 = Carbon('CG2')
        self.add_atom(self.CG2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.OG1))
        self.add_bond(Bond(self.CB, self.CG2))
