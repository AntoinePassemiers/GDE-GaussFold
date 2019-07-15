# -*- coding: utf-8 -*-
# glutamine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Group, Carbon, Hydrogen, Oxygen, Nitrogen


class Glutamine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'GLN', 'Q', **kwargs)

        # Add side chain atoms

        self.COD_group = Group('COD')
        self.NE2_group = Group('NE2')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD = Carbon('CD', q=0.38)
        self.add_atom(self.CD)
        self.COD_group.add(self.CD)

        self.OE1 = Oxygen('OE1', q=-0.38)
        self.add_atom(self.OE1)
        self.COD_group.add(self.OE1)

        self.NE2 = Nitrogen('NE2', q=-0.56)
        self.add_atom(self.NE2)
        self.NE2_group.add(self.NE2)

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

        self.HNA = Hydrogen('HNA', q=0.28)
        self.add_atom(self.HNA)
        self.add_bond(Bond(self.HNA, self.NE2))
        self.NE2_group.add(self.HNA)

        self.HNB = Hydrogen('HNB', q=0.28)
        self.add_atom(self.HNB)
        self.add_bond(Bond(self.HNB, self.NE2))
        self.NE2_group.add(self.HNB)
