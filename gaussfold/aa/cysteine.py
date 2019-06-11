# -*- coding: utf-8 -*-
# cysteine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Sulfur


class Cysteine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'CYS', 'C')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.SG = Sulfur('SG')
        self.add_atom(self.SG)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.SG, self.CB))

        # TODO: allow disulfide bond between paired cysteines
