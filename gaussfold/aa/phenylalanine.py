# -*- coding: utf-8 -*-
# phenylalanine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen


class Phenylalanine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'PHE', 'F')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD1 = Carbon('CD1')
        self.add_atom(self.CD1)

        self.CD2 = Carbon('CD2')
        self.add_atom(self.CD2)

        self.CE1 = Carbon('CE1')
        self.add_atom(self.CE1)

        self.CE2 = Carbon('CE2')
        self.add_atom(self.CE2)

        self.CZ = Carbon('CZ')
        self.add_atom(self.CZ)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.CD1, self.CG, order=2))
        self.add_bond(Bond(self.CD2, self.CG))
        self.add_bond(Bond(self.CE1, self.CD1))
        self.add_bond(Bond(self.CE2, self.CD2, order=2))
        self.add_bond(Bond(self.CZ, self.CE1, order=2))
