# -*- coding: utf-8 -*-
# carbon.py: Carbon atoms
# author : Antoine Passemiers

from gaussfold.atom.atom import Atom


class Carbon(Atom):

    def __init__(self, name, **kwargs):
        Atom.__init__(self, 'C', name, **kwargs)
