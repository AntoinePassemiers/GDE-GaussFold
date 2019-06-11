# -*- coding: utf-8 -*-
# methionine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen, Sulfur


class Methionine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'MET', 'M')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.SD = Sulfur('SD')
        self.add_atom(self.SD)

        self.CE = Carbon('CE')
        self.add_atom(self.CE)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.SD, self.CG))
        self.add_bond(Bond(self.CE, self.SD))
