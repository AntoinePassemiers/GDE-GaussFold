# -*- coding: utf-8 -*-
# aspartic_acid.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Group, Carbon, Hydrogen, Oxygen


class AsparticAcid(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'ASP', 'D', **kwargs)

        # Add side chain atoms

        self.ASP_group = Group('ASP')

        self.CB = Carbon('CB', q=-0.20)
        self.add_atom(self.CB)
        self.ASP_group.add(self.CB)

        self.CG = Carbon('CG', q=0.34)
        self.add_atom(self.CG)
        self.ASP_group.add(self.CG)

        self.OD1 = Oxygen('OD1', q=-0.57)
        self.add_atom(self.OD1)
        self.ASP_group.add(self.OD1)

        self.OD2 = Oxygen('OD2', q=-0.57)
        self.add_atom(self.OD2)
        self.ASP_group.add(self.OD2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.OD1, self.CG, order=2)) # TODO
        self.add_bond(Bond(self.OD2, self.CG)) # TODO

        # Add hydrogens

        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CB))

        self.H2 = Hydrogen('H2')
        self.add_atom(self.H2)
        self.add_bond(Bond(self.H2, self.CB))

        self.H3 = Hydrogen('H3')
        self.add_atom(self.H3)
        self.add_bond(Bond(self.H3, self.OD2))
