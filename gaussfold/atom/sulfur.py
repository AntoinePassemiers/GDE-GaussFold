# -*- coding: utf-8 -*-
# sulfur.py: Sulfur atoms
# author : Antoine Passemiers

from gaussfold.atom.atom import Atom


class Sulfur(Atom):

    def __init__(self, name, **kwargs):
        Atom.__init__(self, 'S', name, **kwargs)
