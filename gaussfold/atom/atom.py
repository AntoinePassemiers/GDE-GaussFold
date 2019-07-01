# -*- coding: utf-8 -*-
# atom.py: Base class for atoms
# author : Antoine Passemiers

import numpy as np
from abc import ABCMeta, abstractmethod


class Atom(metaclass=ABCMeta):

    def __init__(self, element, name, coords=(0., 0., 0.), identifier=None):
        self.element = element
        self.name = name
        self.identifier = identifier
        self.bonded_atoms = list()
        self.coords = np.asarray(list(coords), dtype=np.float)

    def add_bonded_atom(self, atom):
        assert(isinstance(atom, Atom))
        if atom not in self.bonded_atoms:
            self.bonded_atoms.append(atom)

    def set_coords(self, x, y, z):
        self.coords[0] = x
        self.coords[1] = y
        self.coords[2] = z

    def get_coords(self):
        return self.coords

    def is_bonded(self, atom):
        return (atom in self.bonded_atoms)

    def __to_pdb__(self, serial, chain_id, res_seq):
        s = 'ATOM  '
        s += str(serial).rjust(5) + '  '
        s += str(self.name).ljust(3) + ' '
        s += chain_id[0]
        s += str(res_seq).rjust(4) + '    '
        s += str(self.x).ljust(8)
        s += str(self.y).ljust(8)
        s += str(self.z).ljust(8)
        s += '1.00  ' # TODO
        s += '1.00' # TODO
        s += '           '
        s += self.element
        s += '\n'
        return s

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]

    @property
    def z(self):
        return self.coords[2]
