# -*- coding: utf-8 -*-
# atom.py: Base class for atoms
# author : Antoine Passemiers

from gaussfold.atom.group import Group

import numpy as np
from abc import ABCMeta, abstractmethod


class Atom(metaclass=ABCMeta):

    __n_atoms__ = 0

    def __init__(self, element, name, coords=(np.nan, np.nan, np.nan), q=0.):
        self.element = element
        self.name = name
        self._identifier = Atom.__n_atoms__
        Atom.__n_atoms__ += 1
        self.bonded_atoms = list()
        self.coords = np.asarray(list(coords), dtype=np.float)
        self.charge = q
        
        self._group = Group(self.name)
        self._group.add(self)

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

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]

    @property
    def z(self):
        return self.coords[2]

    def __hash__(self):
        return self._identifier

    def __to_pdb__(self, serial, chain_id, res_seq):
        s = 'ATOM  '
        s += str(serial).rjust(5) + '  '
        s += str(self.name).ljust(3) + ' '
        s += '    ' # TODO: res name
        s += chain_id[0] + ' '
        s += str(res_seq).ljust(4) + '   '
        s += '{0:.3f}'.format(self.x).rjust(8)
        s += '{0:.3f}'.format(self.y).rjust(8)
        s += '{0:.3f}'.format(self.z).rjust(8)
        s += '  1.00  ' # TODO
        s += '1.00' # TODO
        s += '           '
        s += self.element
        return s
