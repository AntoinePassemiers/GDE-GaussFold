# -*- coding: utf-8 -*-
# gaussian_constraint.py
# author : Antoine Passemiers

from gaussfold.atom.dummy import DummyAtom

from abc import abstractmethod, ABCMeta


class GaussianConstraint(metaclass=ABCMeta):

    __CENTER_OF_MASS__ = DummyAtom('center of mass')

    def __init__(self, atom_a, atom_b, weight=1.):
        self.atom_a = atom_a
        self.atom_b = atom_b
        self._weight = weight

    def atoms(self):
        return (self.atom_a, self.atom_b)

    @abstractmethod
    def mu(self):
        pass

    @abstractmethod
    def sigma(self):
        pass

    def weight(self):
        return self._weight
