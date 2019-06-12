# -*- coding: utf-8 -*-
# example.py: GDE-GaussFold example
# author : Antoine Passemiers

from gaussfold import GaussFold, Optimizer, tm_score
from gaussfold import PDBParser, SS3Parser, FastaParser, ContactParser

import os
import numpy as np


DATA_FOLDER = 'data'


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

    # Parse native 3D structure from PDB file
    # Residue coordinates are based on C-alpha atoms
    sequence_name = '1DBRA'
    filepath = os.path.join(DATA_FOLDER, 'native.pdb')
    parser = PDBParser(sequence, sequence_name, method='CA')
    distances, coords_target = parser.parse(filepath)

    # Parse predicted contact probabilities
    L = len(sequence) # Number of residues
    filepath = os.path.join(DATA_FOLDER, 'predicted.con')
    #cmap = ContactParser(L).parse(filepath)
    cmap = 1. / (1 + np.exp(distances - 8.))

    # Create GDE-GaussFold object
    gf = GaussFold(n_init_sols=1)

    # Set optimizer hyper-parameters (optional)
    gf.optimizer = Optimizer(
        pop_size=2000,       # Population size
        n_iter=1000,        # Maximum number of iterations
        partition_size=50,   # Partition size for the selection of parents
        mutation_rate=0.5,   # Percentage of child's points to be mutated
        mutation_std=0.3,    # Stdv of mutation noise
        init_std=10.,        # Stdv for randomly generating initial solutions
        early_stopping=5000) # Maximum number of iterations without improvement

    # Run GDE-GaussFold
    coords_predicted = gf.run(cmap, ssp, sequence, verbose=True)

    # Make 3D alignement between native structure and predicted structure
    # and compute TM-score
    print('Compute TM-score')
    tm, coords_predicted = tm_score(coords_predicted, coords_target, return_coords=True)
    print('tm: ', tm)

    import numpy as np
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    coords_target = np.asarray([x for x in coords_target if x is not None])
    x = coords_target[:, 0]
    y = coords_target[:, 1]
    z = coords_target[:, 2]
    ax.plot(x, y, z)
    x = coords_predicted[:, 0]
    y = coords_predicted[:, 1]
    z = coords_predicted[:, 2]
    ax.plot(x, y, z)
    ax.legend()
    plt.show()
