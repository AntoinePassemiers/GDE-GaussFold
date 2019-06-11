# -*- coding: utf-8 -*-
# lysine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Lysine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'LYS', 'K')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD = Carbon('CD')
        self.add_atom(self.CD)

        self.CE = Carbon('CE')
        self.add_atom(self.CE)

        self.NZ = Nitrogen('NZ')
        self.add_atom(self.NZ)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.CD, self.CG))
        self.add_bond(Bond(self.CE, self.CD))
        self.add_bond(Bond(self.NZ, self.CE))
