# -*- coding: utf-8 -*-
# oxygen.py: Oxygen atoms
# author : Antoine Passemiers

from gaussfold.atom.atom import Atom


class Oxygen(Atom):

    def __init__(self, name, **kwargs):
        Atom.__init__(self, 'O', name, **kwargs)
