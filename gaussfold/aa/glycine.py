# -*- coding: utf-8 -*-
# glycine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Hydrogen


class Glycine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'GLY', 'G', **kwargs)

        # Add hydrogen
        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CA))
