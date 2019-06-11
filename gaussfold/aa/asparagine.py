# -*- coding: utf-8 -*-
# asparagine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Asparagine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'ASN', 'N')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.OD1 = Oxygen('OD1')
        self.add_atom(self.OD1)

        self.ND2 = Nitrogen('ND2')
        self.add_atom(self.ND2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.OD1, self.CG, order=2))
        self.add_bond(Bond(self.ND2, self.CG))
