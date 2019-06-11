# -*- coding: utf-8 -*-
# valine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Valine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'VAL')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG1 = Carbon('CG1')
        self.add_atom(self.CG1)

        self.CG2 = Carbon('CG2')
        self.add_atom(self.CG2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.CG1))
        self.add_bond(Bond(self.CB, self.CG2))
