# -*- coding: utf-8 -*-
# hydrogen.py: Hydrogen atoms
# author : Antoine Passemiers

from gaussfold.atom.atom import Atom


class Hydrogen(Atom):

    def __init__(self, name, **kwargs):
        Hydrogen.__init__(self, 'H', name, **kwargs)
