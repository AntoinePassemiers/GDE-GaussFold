# -*- coding: utf-8 -*-
# nitrogen.py: Nitrogen atoms
# author : Antoine Passemiers

from gaussfold.atom.atom import Atom


class Nitrogen(Atom):

    def __init__(self, name, **kwargs):
        Atom.__init__(self, 'N', name, **kwargs)
