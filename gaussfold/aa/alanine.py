# -*- coding: utf-8 -*-
# alaline.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Alanine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'ALA')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.add_bond(Bond(self.CB, self.CA))
