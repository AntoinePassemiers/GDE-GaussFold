# -*- coding: utf-8 -*-
# threonine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Hydrogen, Oxygen, Nitrogen


class Threonine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'THR', 'T', **kwargs)

        # Add side chain atoms

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.OG1 = Oxygen('OG1')
        self.add_atom(self.OG1)

        self.CG2 = Carbon('CG2')
        self.add_atom(self.CG2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.OG1))
        self.add_bond(Bond(self.CB, self.CG2))

        # Add hydrogens

        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CB))

        self.H2 = Hydrogen('H2')
        self.add_atom(self.H2)
        self.add_bond(Bond(self.H2, self.OG1))

        self.H3 = Hydrogen('H3')
        self.add_atom(self.H3)
        self.add_bond(Bond(self.H3, self.CG2))

        self.H4 = Hydrogen('H4')
        self.add_atom(self.H4)
        self.add_bond(Bond(self.H4, self.CG2))

        self.H5 = Hydrogen('H5')
        self.add_atom(self.H5)
        self.add_bond(Bond(self.H5, self.CG2))
