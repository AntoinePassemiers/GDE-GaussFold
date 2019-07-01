# -*- coding: utf-8 -*-
# atom.py: Base class for atoms
# author : Antoine Passemiers

from abc import ABCMeta, abstractmethod


class Atom(metaclass=ABCMeta):

    def __init__(self, element, name, coords=(0., 0., 0.), identifier=None):
        self.element = element
        self.name = name
        self.identifier = identifier
        self.bonded_atoms = list()

        self.x = float(coords[0])
        self.y = float(coords[1])
        self.z = float(coords[2])

    def add_bonded_atom(self, atom):
        assert(isinstance(atom, Atom))
        if atom not in self.bonded_atoms:
            self.bonded_atoms.append(atom)

    def set_coords(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def get_coords(self):
        return (self.x, self.y, self.z)

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
