# -*- coding: utf-8 -*-
# arginine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Arginine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'ARG')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD = Carbon('CD')
        self.add_atom(self.CD)

        self.NE = Nitrogen('NE')
        self.add_atom(self.NE)

        self.CZ = Carbon('CZ')
        self.add_atom(self.CZ)

        self.NH1 = Nitrogen('NH1')
        self.add_atom(self.NH1)

        self.NH2 = Nitrogen('NH2')
        self.add_atom(self.NH2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.CG))
        self.add_bond(Bond(self.CG, self.CD))
        self.add_bond(Bond(self.CD, self.NE))
        self.add_bond(Bond(self.NE, self.CZ))
        self.add_bond(Bond(self.CZ, self.NH1))
        self.add_bond(Bond(self.CZ, self.NH2))
