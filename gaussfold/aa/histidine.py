# -*- coding: utf-8 -*-
# histidine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Histidine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'HIS')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.ND1 = Nitrogen('ND1')
        self.add_atom(self.ND1)

        self.CD2 = Carbon('CD2')
        self.add_atom(self.CD2)

        self.CE1 = Carbon('CE1')
        self.add_atom(self.CE1)

        self.NE2 = Nitrogen('NE2')
        self.add_atom(self.NE2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.CG))
        self.add_bond(Bond(self.CG, self.ND1))
        self.add_bond(Bond(self.CG, self.CD2, order=2))
        self.add_bond(Bond(self.CD2, self.NE2))
        self.add_bond(Bond(self.ND1, self.CE1, order=2))
        self.add_bond(Bond(self.NE2, self.CE1))
