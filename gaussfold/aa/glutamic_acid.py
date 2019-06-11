# -*- coding: utf-8 -*-
# glutamic_acid.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class GlutamicAcid(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'GLU', 'E')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD = Carbon('CD')
        self.add_atom(self.CD)

        self.OE1 = Oxygen('OE1')
        self.add_atom(self.OE1)

        self.OE2 = Oxygen('OE2')
        self.add_atom(self.OE2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.CD, self.CG))
        self.add_bond(Bond(self.OE1, self.CD, order=2)) # TODO
        self.add_bond(Bond(self.OE2, self.CD)) # TODO
