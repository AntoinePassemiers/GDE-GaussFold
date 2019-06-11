# -*- coding: utf-8 -*-
# glutamine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Glutamine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'GLN', 'Q')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD = Carbon('CD')
        self.add_atom(self.CD)

        self.OE1 = Oxygen('OE1')
        self.add_atom(self.OE1)

        self.NE2 = Nitrogen('NE2')
        self.add_atom(self.NE2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.CD, self.CG))
        self.add_bond(Bond(self.OE1, self.CD, order=2))
        self.add_bond(Bond(self.NE2, self.CD))
