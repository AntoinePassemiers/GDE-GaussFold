# -*- coding: utf-8 -*-
# histidine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Group, Carbon, Hydrogen, Nitrogen


class Histidine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'HIS', 'H', **kwargs)

        # Add side chain atoms

        self.HIS_group = Group('HIS')

        self.CB = Carbon('CB', q=0.11)
        self.add_atom(self.CB)
        self.HIS_group.add(self.CB)

        self.CG = Carbon('CG', q=0.17)
        self.add_atom(self.CG)
        self.HIS_group.add(self.CG)

        self.ND1 = Nitrogen('ND1', q=-0.50)
        self.add_atom(self.ND1)
        self.HIS_group.add(self.ND1)

        self.CD2 = Carbon('CD2', q=0.33)
        self.add_atom(self.CD2)
        self.HIS_group.add(self.CD2)

        self.CE1 = Carbon('CE1', q=0.65)
        self.add_atom(self.CE1)
        self.HIS_group.add(self.CE1)

        self.NE2 = Nitrogen('NE2', q=-0.50)
        self.add_atom(self.NE2)
        self.HIS_group.add(self.NE2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.CG))
        self.add_bond(Bond(self.CG, self.ND1))
        self.add_bond(Bond(self.CG, self.CD2, order=2))
        self.add_bond(Bond(self.CD2, self.NE2))
        self.add_bond(Bond(self.ND1, self.CE1, order=2))
        self.add_bond(Bond(self.NE2, self.CE1))

        # Add hydrogens

        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CB))

        self.H2 = Hydrogen('H2')
        self.add_atom(self.H2)
        self.add_bond(Bond(self.H2, self.CB))

        self.H3 = Hydrogen('H3')
        self.add_atom(self.H3)
        self.add_bond(Bond(self.H3, self.CD2))

        self.HE2 = Hydrogen('HE2', q=0.37)
        self.add_atom(self.HE2)
        self.add_bond(Bond(self.HE2, self.NE2))
        self.HIS_group.add(self.HE2)

        self.H5 = Hydrogen('H5')
        self.add_atom(self.H5)
        self.add_bond(Bond(self.H5, self.CE1))
