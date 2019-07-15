# -*- coding: utf-8 -*-
# group.py: Atom groups
# author : Antoine Passemiers


class Group:

    def __init__(self, name):
        self._name = name
        self._atoms = list()

    def add(self, atom):
        atom._group = self
        self._atoms.append(atom)
