# -*- coding: utf-8 -*-
# glycine.py
# author : Antoine Passemiers

from gaussfold.aa.amino_acid import AminoAcid
from gaussfold.atom import Bond, Carbon, Oxygen, Nitrogen


class Glycine(AminoAcid):

    def __init__(self):
        AminoAcid.__init__(self, 'GLY')
