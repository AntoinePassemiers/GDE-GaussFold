# -*- coding: utf-8 -*-
# glutamine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Hydrogen, Oxygen, Nitrogen


class Glutamine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'GLN', 'Q', **kwargs)

        # Add side chain atoms

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD = Carbon('CD')
        self.add_atom(self.CD)

        self.OE1 = Oxygen('OE1')
        self.add_atom(self.OE1)

        self.NE2 = Nitrogen('NE2')
        self.add_atom(self.NE2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.CD, self.CG))
        self.add_bond(Bond(self.OE1, self.CD, order=2))
        self.add_bond(Bond(self.NE2, self.CD))

        # Add hydrogens

        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CB))

        self.H2 = Hydrogen('H2')
        self.add_atom(self.H2)
        self.add_bond(Bond(self.H2, self.CB))

        self.H3 = Hydrogen('H3')
        self.add_atom(self.H3)
        self.add_bond(Bond(self.H3, self.CG))

        self.H4 = Hydrogen('H4')
        self.add_atom(self.H4)
        self.add_bond(Bond(self.H4, self.CG))

        self.H5 = Hydrogen('H5')
        self.add_atom(self.H5)
        self.add_bond(Bond(self.H5, self.NE2))

        self.H6 = Hydrogen('H6')
        self.add_atom(self.H6)
        self.add_bond(Bond(self.H6, self.NE2))
