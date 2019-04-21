# -*- coding: utf-8 -*-
# optimizer.py
# author : Antoine Passemiers

import random
import numpy as np


class Optimizer:

    def __init__(self, initial_coords, obj):
        self.L = len(initial_coords)
        self.initial_coords = initial_coords
        self.obj = obj
        self.scores = list()

    def random_sol(self):
        offsets = np.random.normal(0., 10., size=(self.L, 3))
        sol = self.initial_coords + offsets
        return sol

    def cross_over(self, left, right):
        alpha = np.random.randint(0, 2, size=(len(left), 1))
        s = alpha * left + (1. - alpha) * right
        mutations = np.random.normal(0., 0.3, size=(self.L, 3))
        mutations *= np.random.randint(0, 2, size=(len(left), 1))
        individual = s + mutations
        return individual

    def new_sol(self, pop, scores, ps=50):
        pop = [_ for _ in pop]
        np.copy(scores)
        indices = np.arange(len(pop))
        random.shuffle(indices)
        pop = [pop[i] for i in indices]
        scores = scores[indices]

        left_winner = pop[np.argmax(scores[:ps])]
        right_winner = pop[ps + np.argmax(scores[ps:])]

        return self.cross_over(left_winner, right_winner)

    def run(self, pop_size=2000, n_iter=100000, partition_size=50,
            early_stopping=2, verbose=True):
        pop = [self.random_sol() for i in range(pop_size-1)]
        pop.append(self.initial_coords)

        self.scores = list()
        best_score = self.obj(pop[-1])

        best_iteration = 0

        scores = np.asarray([self.obj(ind) for ind in pop])
        for k in range(n_iter):
            new_ind = self.new_sol(pop, scores, ps=partition_size)
            worst = np.argmin(scores)
            pop[worst] = new_ind
            scores[worst] = self.obj(new_ind)

            if scores[worst] > best_score:
                best_score = scores[worst]
                best_iteration = k
            if verbose and (k + 1) % 100 == 0:
                print('Log-likelihood at iteration %i: %f' \
                    % (k + 1, best_score))
            self.scores.append(best_score)

            if k - best_iteration >= early_stopping:
                break

        best_coords = pop[np.argmax(scores)]
        return best_coords
