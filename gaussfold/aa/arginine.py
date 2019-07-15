# -*- coding: utf-8 -*-
# arginine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Group, Carbon, Hydrogen, Nitrogen


class Arginine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'ARG', 'R', **kwargs)

        # Add side chain atoms

        self.ARG_group = Group('ARG')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD = Carbon('CD', q=0.19)
        self.add_atom(self.CD)
        self.ARG_group.add(self.CD)

        self.NE = Nitrogen('NE', q=-0.50)
        self.add_atom(self.NE)
        self.ARG_group.add(self.NE)

        self.CZ = Carbon('CZ', q=0.46)
        self.add_atom(self.CZ)
        self.ARG_group.add(self.CZ)

        self.NH1 = Nitrogen('NH1', q=-0.50)
        self.add_atom(self.NH1)
        self.ARG_group.add(self.NH1)

        self.NH2 = Nitrogen('NH2', q=-0.50)
        self.add_atom(self.NH2)
        self.ARG_group.add(self.NH2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CB, self.CG))
        self.add_bond(Bond(self.CG, self.CD))
        self.add_bond(Bond(self.CD, self.NE))
        self.add_bond(Bond(self.NE, self.CZ))
        self.add_bond(Bond(self.CZ, self.NH1, order=2)) # TODO
        self.add_bond(Bond(self.CZ, self.NH2)) # TODO

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
        self.add_bond(Bond(self.H5, self.CD))

        self.H6 = Hydrogen('H6')
        self.add_atom(self.H6)
        self.add_bond(Bond(self.H6, self.CD))

        self.HNE = Hydrogen('HNE', q=0.37)
        self.add_atom(self.HNE)
        self.add_bond(Bond(self.HNE, self.NE))
        self.ARG_group.add(self.HNE)

        self.HHA = Hydrogen('HHA', q=0.37)
        self.add_atom(self.HHA)
        self.add_bond(Bond(self.HHA, self.NH1))
        self.ARG_group.add(self.HHA)

        self.HHB = Hydrogen('HHB', q=0.37)
        self.add_atom(self.HHB)
        self.add_bond(Bond(self.HHB, self.NH2))
        self.ARG_group.add(self.HHB)

        self.HHC = Hydrogen('HHC', q=0.37)
        self.add_atom(self.HHC)
        self.add_bond(Bond(self.HHC, self.NH2))
        self.ARG_group.add(self.HHC)
