# -*- coding: utf-8 -*-
# serine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Group, Carbon, Hydrogen, Oxygen, Nitrogen


class Serine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'SER', 'S', **kwargs)

        # Add side chain atoms

        self.OGS_group = Group('OGS')

        self.CB = Carbon('CB', q=0.03)
        self.add_atom(self.CB)
        self.OGS_group.add(self.CB)

        self.OG = Oxygen('OG', q=-0.38)
        self.add_atom(self.OG)
        self.OGS_group.add(self.OG)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.OG))

        # Add hydrogens

        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CB))

        self.H2 = Hydrogen('H2')
        self.add_atom(self.H2)
        self.add_bond(Bond(self.H2, self.CB))

        self.HOG = Hydrogen('HOG', q=-0.35)
        self.add_atom(self.HOG)
        self.add_bond(Bond(self.HOG, self.OG))
        self.OGS_group.add(self.HOG)
