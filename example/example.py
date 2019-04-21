# -*- coding: utf-8 -*-
# example.py: GDE-GaussFold example
# author : Antoine Passemiers

from gaussfold import GaussFold, tm_score, PDBParser, SS3Parser, FastaParser, ContactParser

import os


DATA_FOLDER = 'data'


if __name__ == '__main__':

    sequence = FastaParser().parse(
        os.path.join(DATA_FOLDER, 'sequence.fa'))['sequences'][0]
    print(sequence)

    ssp = SS3Parser().parse(
        os.path.join(DATA_FOLDER, 'ss3.txt')).argmax(axis=1)
    print(ssp)

    sequence_name = '1DBRA'
    distances, coords_target = PDBParser(sequence, sequence_name,
        method='CA').parse(os.path.join(DATA_FOLDER, 'native.pdb'))

    L = len(ssp)
    cmap = ContactParser(L).parse(
        os.path.join(DATA_FOLDER, 'predicted.con'))

    coords_predicted = GaussFold().run(cmap, ssp)

    print('Compute TM-score')
    print('tm: ', tm_score(coords_predicted, coords_target))
