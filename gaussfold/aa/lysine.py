# -*- coding: utf-8 -*-
# lysine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Group, Carbon, Hydrogen, Nitrogen


class Lysine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'LYS', 'K', **kwargs)

        # Add side chain atoms

        self.LYS_group = Group('LYS')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD = Carbon('CD', q=0.12)
        self.add_atom(self.CD)
        self.LYS_group.add(self.CD)

        self.CE = Carbon('CE', q=0.30)
        self.add_atom(self.CE)
        self.LYS_group.add(self.CE)

        self.NZ = Nitrogen('NZ', q=-0.50)
        self.add_atom(self.NZ)
        self.LYS_group.add(self.NZ)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.CD, self.CG))
        self.add_bond(Bond(self.CE, self.CD))
        self.add_bond(Bond(self.NZ, self.CE))

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

        self.H7 = Hydrogen('H7')
        self.add_atom(self.H7)
        self.add_bond(Bond(self.H7, self.CE))

        self.H8 = Hydrogen('H8')
        self.add_atom(self.H8)
        self.add_bond(Bond(self.H8, self.CE))

        self.HZA = Hydrogen('HZA', q=0.36)
        self.add_atom(self.HZA)
        self.add_bond(Bond(self.HZA, self.NZ))
        self.LYS_group.add(self.HZA)

        self.HZB = Hydrogen('HZB', q=0.36)
        self.add_atom(self.HZB)
        self.add_bond(Bond(self.HZB, self.NZ))
        self.LYS_group.add(self.HZB)
