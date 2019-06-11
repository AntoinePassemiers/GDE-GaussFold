# -*- coding: utf-8 -*-
# serine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Serine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'SER')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.OG = Oxygen('OG')
        self.add_atom(self.OG)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.OG))
