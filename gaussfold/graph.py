# -*- coding: utf-8 -*-
# graph.py: GDE-GaussFold core algorithm
# author : Antoine Passemiers

import numpy as np
import networkx as nx
from networkx.algorithms.shortest_paths.generic import shortest_path_length


class Graph:

    def __init__(self, A):
        self.A = np.asarray(A)
        self.L = self.A.shape[0]
        self.G = nx.from_numpy_array(self.A, parallel_edges=False)

    def is_connected(self):
        return nx.is_connected(self.G)

    def distances(self):
        path_lengths = shortest_path_length(self.G)
        gds = np.zeros((self.L, self.L), dtype=np.int)
        for i, i_lengths in path_lengths:
            for j in i_lengths.keys():
                gds[i, j] = gds[j, i] = i_lengths[j]
        return gds
