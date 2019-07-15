# -*- coding: utf-8 -*-
# tyrosine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Group, Carbon, Hydrogen, Oxygen


class Tyrosine(AminoAcid):

    def __init__(self, **kwargs):
        AminoAcid.__init__(self, 'TYR', 'Y', **kwargs)

        # Add side chain atoms

        self.OHH_group = Group('OHH')

        self.CB = Carbon('CB')
        self.add_atom(self.CB)

        self.CG = Carbon('CG')
        self.add_atom(self.CG)

        self.CD1 = Carbon('CD1')
        self.add_atom(self.CD1)

        self.CD2 = Carbon('CD2')
        self.add_atom(self.CD2)

        self.CE1 = Carbon('CE1')
        self.add_atom(self.CE1)

        self.CE2 = Carbon('CE2')
        self.add_atom(self.CE2)

        self.CZ = Carbon('CZ', q=0.03)
        self.add_atom(self.CZ)
        self.OHH_group.add(self.CZ)

        self.OH = Oxygen('OH', q=-0.38)
        self.add_atom(self.OH)
        self.OHH_group.add(self.OH)

        self.add_bond(Bond(self.CB, self.CA))
        self.add_bond(Bond(self.CG, self.CB))
        self.add_bond(Bond(self.CD1, self.CG, order=2))
        self.add_bond(Bond(self.CD2, self.CG))
        self.add_bond(Bond(self.CE1, self.CD1))
        self.add_bond(Bond(self.CE2, self.CD2, order=2))
        self.add_bond(Bond(self.CZ, self.CE1, order=2))
        self.add_bond(Bond(self.OH, self.CZ))

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
        self.add_bond(Bond(self.H4, self.CD2))

        self.H5 = Hydrogen('H5')
        self.add_atom(self.H5)
        self.add_bond(Bond(self.H5, self.CE1))

        self.H6 = Hydrogen('H6')
        self.add_atom(self.H6)
        self.add_bond(Bond(self.H6, self.CE2))

        self.HOH = Hydrogen('HOH', q=0.35)
        self.add_atom(self.HOH)
        self.add_bond(Bond(self.HOH, self.OH))
        self.OHH_group.add(self.HOH)
