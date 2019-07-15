# -*- coding: utf-8 -*-
# threonine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Group, Carbon, Hydrogen, Oxygen, Nitrogen


class Threonine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'THR', 'T', **kwargs)

        # Add side chain atoms

        self.OGT_group = Group('OGT')

        self.CB = Carbon('CB', q=0.03)
        self.add_atom(self.CB)
        self.OGT_group.add(self.CB)

        self.OG1 = Oxygen('OG1', q=-0.38)
        self.add_atom(self.OG1)
        self.OGT_group.add(self.OG1)

        self.CG2 = Carbon('CG2')
        self.add_atom(self.CG2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.OG1))
        self.add_bond(Bond(self.CB, self.CG2))

        # Add hydrogens

        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CB))

        self.HOG = Hydrogen('HOG', q=0.35)
        self.add_atom(self.HOG)
        self.add_bond(Bond(self.HOG, self.OG1))
        self.OGT_group.add(self.HOG)

        self.H3 = Hydrogen('H3')
        self.add_atom(self.H3)
        self.add_bond(Bond(self.H3, self.CG2))

        self.H4 = Hydrogen('H4')
        self.add_atom(self.H4)
        self.add_bond(Bond(self.H4, self.CG2))

        self.H5 = Hydrogen('H5')
        self.add_atom(self.H5)
        self.add_bond(Bond(self.H5, self.CG2))
