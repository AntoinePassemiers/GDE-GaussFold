# -*- coding: utf-8 -*-
# arginine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Hydrogen, Nitrogen


class Arginine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'ARG', 'R', **kwargs)

        # Add side chain atoms

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD = Carbon('CD')
        self.add_atom(self.CD)

        self.NE = Nitrogen('NE')
        self.add_atom(self.NE)

        self.CZ = Carbon('CZ')
        self.add_atom(self.CZ)

        self.NH1 = Nitrogen('NH1')
        self.add_atom(self.NH1)

        self.NH2 = Nitrogen('NH2')
        self.add_atom(self.NH2)

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

        self.H7 = Hydrogen('H7')
        self.add_atom(self.H7)
        self.add_bond(Bond(self.H7, self.NE))

        self.H8 = Hydrogen('H8')
        self.add_atom(self.H8)
        self.add_bond(Bond(self.H8, self.NH1))

        self.H9 = Hydrogen('H9')
        self.add_atom(self.H9)
        self.add_bond(Bond(self.H9, self.NH2))

        self.H10 = Hydrogen('H10')
        self.add_atom(self.H10)
        self.add_bond(Bond(self.H10, self.NH2))
