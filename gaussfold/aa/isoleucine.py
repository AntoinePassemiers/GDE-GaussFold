# -*- coding: utf-8 -*-
# isoleucine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Isoleucine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'ILE', 'I')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG1 = Carbon('CG1')
        self.add_atom(self.CG1)

        self.CG2 = Carbon('CG2')
        self.add_atom(self.CG2)

        self.CD1 = Carbon('CD1')
        self.add_atom(self.CD1)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG1, self.CB))
        self.add_bond(Bond(self.CG2, self.CB))
        self.add_bond(Bond(self.CD1, self.CG1))
