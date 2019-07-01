# -*- coding: utf-8 -*-
# leucine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Hydrogen


class Leucine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'LEU', 'L', **kwargs)

        # Add side chain atoms

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD1 = Carbon('CD1')
        self.add_atom(self.CD1)

        self.CD2 = Carbon('CD2')
        self.add_atom(self.CD2)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.CD1, self.CG))
        self.add_bond(Bond(self.CD2, self.CG))

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
        self.add_bond(Bond(self.H4, self.CD2))

        self.H5 = Hydrogen('H5')
        self.add_atom(self.H5)
        self.add_bond(Bond(self.H5, self.CD2))

        self.H6 = Hydrogen('H6')
        self.add_atom(self.H6)
        self.add_bond(Bond(self.H6, self.CD2))

        self.H7 = Hydrogen('H7')
        self.add_atom(self.H7)
        self.add_bond(Bond(self.H7, self.CD1))

        self.H8 = Hydrogen('H8')
        self.add_atom(self.H8)
        self.add_bond(Bond(self.H8, self.CD1))

        self.H9 = Hydrogen('H9')
        self.add_atom(self.H9)
        self.add_bond(Bond(self.H9, self.CD1))
