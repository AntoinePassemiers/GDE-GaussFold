# -*- coding: utf-8 -*-
# tryptophan.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Hydrogen, Nitrogen


class Tryptophan(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'TRP', 'W', **kwargs)

        # Add side chain atoms

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD1 = Carbon('CD1')
        self.add_atom(self.CD1)

        self.CD2 = Carbon('CD2')
        self.add_atom(self.CD2)

        self.NE1 = Nitrogen('NE1')
        self.add_atom(self.NE1)

        self.CE2 = Carbon('CE2')
        self.add_atom(self.CE2)

        self.CE3 = Carbon('CE3')
        self.add_atom(self.CE3)

        self.CZ2 = Carbon('CZ2')
        self.add_atom(self.CZ2)

        self.CZ3 = Carbon('CZ3')
        self.add_atom(self.CZ3)

        self.CH2 = Carbon('CH2')
        self.add_atom(self.CH2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.CD1, self.CG, order=2))
        self.add_bond(Bond(self.CD2, self.CG))
        self.add_bond(Bond(self.NE1, self.CD1))
        self.add_bond(Bond(self.NE1, self.CE2))
        self.add_bond(Bond(self.CE2, self.CD2, order=2))
        self.add_bond(Bond(self.CE3, self.CD2))
        self.add_bond(Bond(self.CZ3, self.CE3, order=2))
        self.add_bond(Bond(self.CZ2, self.CE2))
        self.add_bond(Bond(self.CH2, self.CZ2, order=2))
        self.add_bond(Bond(self.CZ3, self.CH2))

        # Add hydrogens

        self.H1 = Hydrogen('H1')
        self.add_atom(self.H1)
        self.add_bond(Bond(self.H1, self.CB))

        self.H2 = Hydrogen('H2')
        self.add_atom(self.H2)
        self.add_bond(Bond(self.H2, self.CB))

        self.H3 = Hydrogen('H3')
        self.add_atom(self.H3)
        self.add_bond(Bond(self.H3, self.CD1))

        self.H4 = Hydrogen('H4')
        self.add_atom(self.H4)
        self.add_bond(Bond(self.H4, self.NE1))

        self.H5 = Hydrogen('H5')
        self.add_atom(self.H5)
        self.add_bond(Bond(self.H5, self.CZ2))

        self.H6 = Hydrogen('H6')
        self.add_atom(self.H6)
        self.add_bond(Bond(self.H6, self.CH2))

        self.H7 = Hydrogen('H7')
        self.add_atom(self.H7)
        self.add_bond(Bond(self.H7, self.CZ3))

        self.H8 = Hydrogen('H8')
        self.add_atom(self.H8)
        self.add_bond(Bond(self.H8, self.CE3))
