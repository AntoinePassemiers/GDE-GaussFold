# -*- coding: utf-8 -*-
# asparagine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Group, Carbon, Hydrogen, Oxygen, Nitrogen


class Asparagine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'ASN', 'N', **kwargs)

        # Add side chain atoms

        self.COG_group = Group('COG')
        self.ND2_group = Group('ND2')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG', q=0.38)
        self.add_atom(self.CG)
        self.COG_group.add(self.CG)

        self.OD1 = Oxygen('OD1', q=-0.38)
        self.add_atom(self.OD1)
        self.COG_group.add(self.OD1)

        self.ND2 = Nitrogen('ND2', q=-0.56)
        self.add_atom(self.ND2)
        self.ND2_group.add(self.ND2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.OD1, self.CG, order=2))
        self.add_bond(Bond(self.ND2, self.CG))

        # Add hydrogens

        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CB))

        self.H2 = Hydrogen('H2')
        self.add_atom(self.H2)
        self.add_bond(Bond(self.H2, self.CB))

        self.HNA = Hydrogen('HNA', q=0.28)
        self.add_atom(self.HNA)
        self.add_bond(Bond(self.HNA, self.ND2))
        self.ND2_group.add(self.HNA)

        self.HNB = Hydrogen('HNB', q=0.28)
        self.add_atom(self.HNB)
        self.add_bond(Bond(self.HNB, self.ND2))
        self.ND2_group.add(self.HNB)
