# -*- coding: utf-8 -*-
# example.py: GDE-GaussFold example
# author : Antoine Passemiers

from gaussfold import GaussFold, Optimizer, tm_score, rmsd
from gaussfold import PDBParser, SS3Parser, FastaParser, ContactParser

import os
import numpy as np


DATA_FOLDER = '../../Wynona/data/additional_set/5DJEB'
#DATA_FOLDER = 'data'


if __name__ == '__main__':

    # Parse primary structure from FASTA file
    filepath = os.path.join(DATA_FOLDER, 'sequence.fa')
    sequence = FastaParser().parse(filepath)['sequences'][0]
    print('\nPrimary structure:')
    print(sequence)

    # Parse predicted 3-state secondary structure
    # and make predictions out of predicted probabilities
    filepath = os.path.join(DATA_FOLDER, 'ss3.txt')
    ssp = SS3Parser().parse(filepath).argmax(axis=1)
    print('\nSecondary structure:')
    print(ssp)
    print('')

    # Parse solvent accesibility probabilities
    filepath = os.path.join(DATA_FOLDER, 'acc.txt')
    acc = SS3Parser().parse(filepath).argmax(axis=1)

    # Parse native 3D structure from PDB file
    # Residue coordinates are based on C-alpha atoms
    # Chain id is assumed to be the last letter of `sequence_name`
    sequence_name = '1DRBA'
    filepath = os.path.join(DATA_FOLDER, 'native.pdb')
    parser = PDBParser(sequence, sequence_name, method='CA')
    distances, coords_target = parser.parse(filepath)

    # Parse predicted contact probabilities
    L = len(sequence) # Number of residues
    filepath = os.path.join(DATA_FOLDER, 'psicov.out')
    if os.path.exists(filepath):
        cmap = ContactParser(L, target_cols=[4]).parse(filepath)
    else:
        filepath = os.path.join(DATA_FOLDER, 'dir.gaussdca')
        cmap = ContactParser(L, target_cols=[2]).parse(filepath)
    #cmap = 1. / (1 + np.exp(distances - 8.))

    # Create GDE-GaussFold object
    gf = GaussFold(n_top=2.5)

    # Set optimizer hyper-parameters (optional)
    gf.optimizer = Optimizer(
        pop_size=2000,       # Population size
        n_iter=50000,        # Maximum number of iterations
        partition_size=50,   # Partition size for the selection of parents
        mutation_rate=0.5,   # Percentage of child's points to be mutated
        mutation_std=0.3,    # Stdv of mutation noise
        init_std=10.,        # Stdv for randomly generating initial solutions
        early_stopping=5000) # Maximum number of iterations without improvement

    # Run GDE-GaussFold
    coords_predicted = gf.run(cmap, ssp, acc, sequence, verbose=True)

    import numpy as np
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import scipy.spatial
    distances_predicted = scipy.spatial.distance.cdist(
            coords_predicted, coords_predicted, metric='euclidean')
    plt.imshow(distances_predicted)
    plt.show()


    # Make 3D alignement between native structure and predicted structure
    # and compute TM-score
    print('Compute TM-score and RMSD')
    tm, coords_aligned = tm_score(coords_predicted, coords_target, return_coords=True)
    print('tm: ', tm)
    r = rmsd(coords_predicted, coords_target, return_coords=False)
    print('r: ', r)


    fig = plt.figure()
    ax = fig.gca(projection='3d')
    coords_target = np.asarray([x for x in coords_target if x is not None])
    x = coords_target[:, 0]
    y = coords_target[:, 1]
    z = coords_target[:, 2]
    ax.plot(x, y, z)
    offset = 25
    x = coords_aligned[:, 0] + offset
    y = coords_aligned[:, 1] + offset
    z = coords_aligned[:, 2] + offset
    ax.plot(x, y, z)
    ax.legend()
    plt.show()
