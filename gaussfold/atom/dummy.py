# -*- coding: utf-8 -*-
# dummy.py
# author : Antoine Passemiers

from gaussfold.atom.atom import Atom


class DummyAtom(Atom):

    def __init__(self, name, **kwargs):
        Atom.__init__(self, '-', name, **kwargs)
