# -*- coding: utf-8 -*-
# chain.py: Base class for amino acid chains
# author : Antoine Passemiers

from gaussfold.aa import *
from gaussfold.aa.arginine import Arginine
from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond

import numpy as np


class Chain:

    def __init__(self, c='CA'):
        self.amino_acids = list()
        assert(c in ['CA', 'CB'])
        self.c = c

    def update(self):
        identifier = 0
        for amino_acid in self.amino_acids:
            for atom in amino_acid.atoms:
                atom.identifier = identifier
                identifier += 1

    def add(self, amino_acid):
        amino_acid.c = self.c
        if len(self.amino_acids) > 0:
            peptide_bond = Bond(self.amino_acids[-1].C, amino_acid.N)
            self.amino_acids[-1].add_bond(peptide_bond)
        self.amino_acids.append(amino_acid)
        self.update()

    def get_atoms_coords(self):
        atoms = self.atoms()
        coords = np.empty((len(atoms), 3), dtype=np.float)
        for i, atom in enumerate(atoms):
            assert(atom.identifier == i)
            coords[i, :] = atom.get_coords()
        return coords

    def set_atoms_coords(self, coords):
        atoms = self.atoms()
        for i, atom in enumerate(atoms):
            assert(atom.identifier == i)
            x, y, z = coords[i, :]
            atom.set_coords(x, y, z)

    def atoms(self):
        atom_list = list()
        for amino_acid in self.amino_acids:
            atom_list += list(amino_acid.atoms)
        return atom_list

    def bonds(self):
        bond_list = list()
        for amino_acid in self.amino_acids:
            bond_list += list(amino_acid.bonds)
        return bond_list

    def dihedral_angles(self):
        phi = np.zeros(self.__len__(), dtype=np.float)
        psi = np.zeros(self.__len__(), dtype=np.float)
        # TODO: terminal amino acids
        for i in range(1, self.__len__() - 1):
            phi[i] = self.amino_acids[i].phi(self.amino_acids[i-1].C)
            psi[i] = self.amino_acids[i].psi(self.amino_acids[i+1].N)
        return phi, psi

    @staticmethod
    def from_string(s, c='CA'):
        assert(c in ('CA', 'CB'))
        abb_to_clk = dict()
        for ckl in AminoAcid.__subclasses__():
            abb_to_clk[ckl(c=c).abbreviation] = ckl

        chain = Chain(c=c)
        for c in s:
            chain.add(abb_to_clk[c]())
        return chain

    def __getitem__(self, key):
        return self.amino_acids[key]

    def __setitem__(self, key, value):
        self.amino_acids[key] = value
        # TODO: check value type
        self.update()

    def __len__(self):
        return len(self.amino_acids)

    def __iter__(self):
        for amino_acid in self.amino_acids:
            yield amino_acid
