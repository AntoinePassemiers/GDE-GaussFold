# -*- coding: utf-8 -*-
# aspartic_acid.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class AsparticAcid(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'ASP', 'D')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.OD1 = Oxygen('OD1')
        self.add_atom(self.OD1)

        self.OD2 = Oxygen('OD2')
        self.add_atom(self.OD2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.OD1, self.CG, order=2)) # TODO
        self.add_bond(Bond(self.OD2, self.CG)) # TODO
