# -*- coding: utf-8 -*-
# optimizer.py: Heuristic optimizer for Gaussian models
# author : Antoine Passemiers

from gaussfold.lbfgs import lbfgs

import random
import numpy as np


class Optimizer:
    """Heuristic optimizer based on a simple vanilla genetic algorithm.

    Attributes:
        pop_size (int): Number of solutions kept in memory.
        n_iter (int): Maximum number of iterations.
        partition_size (int): Partition size for the selection
            of parents. Each partition elects one of the two parents.
        mutation_rate (float): Percentage of points in child's solution
            to be mutated by the given standard deviation.
        mutation_std (float): Standard deviation of the noise to be
            added for mutating points coordinates.
        init_std (float): Standard deviation used to generate the
            population from an initial solution.
        early_stopping (int): Maximum number of iterations without
            score improvement before stopping the algorithm.
        use_lbfgs (bool): Whether to improve local convergence
            with L-BFGS algorithm (slows the solver down).
        scores (list): History of best score over time.
    """

    def __init__(self, pop_size=2000, n_iter=200000, partition_size=50,
                 mutation_rate=0.5, mutation_std=0.3, init_std=10.,
                 early_stopping=300, use_lbfgs=False):
        self.pop_size = pop_size
        self.n_iter = n_iter
        self.partition_size = partition_size
        self.mutation_rate = mutation_rate
        self.mutation_std = mutation_std
        self.init_std = init_std
        self.early_stopping = early_stopping
        self.use_lbfgs = use_lbfgs
        self.scores = list()

    def random_sol(self, initial_coords):
        """Generates a random solution by adding Gaussian noise
        to an initial solution.

        Parameters:
            initial_coords (:obj:`np.ndarray`): Array of shape
                (L, 3) representing the initial solution, where
                L is the number of residues in the protein.

        Returns:
            :obj:`np.ndarray`: A new solution of the same shape.
        """
        L = initial_coords.shape[0]
        offsets = np.random.normal(0., self.init_std, size=(L, 3))
        individual = initial_coords + offsets
        return individual

    def cross_over(self, left, right):
        """Cross-over operator between two parent solutions.

        Parameters:
            left (:obj:`np.ndarray`): First solution of shape (L, 3).
            right (:obj:`np.ndarray`): Second solution of shape (L, 3).

        Returns:
            :obj:`np.ndarray`: Child solution of shape (L, 3).
        """
        alpha = np.random.randint(0, 2, size=(len(left), 1))
        individual = alpha * left + (1. - alpha) * right
        return individual

    def mutate(self, individual):
        """Mutation operator.

        Parameters:
            individual (:obj:`np.ndarray`): Solution of shape (L, 3).

        Returns:
            :obj:`np.ndarray`: Mutated solution of shape (L, 3).
        """
        L = individual.shape[0]
        std = self.mutation_std
        mutations = np.random.normal(0., std, size=(L, 3))
        mutations *= (np.random.rand(L, 1) < self.mutation_rate)
        return individual + mutations

    def new_sol(self, pop, scores):
        """Randomly constructs two partitions from current population,
        plays one tournament in each, elects two parents, applies the
        cross-over operator to get a new solution, and applies the
        mutation operator on it.

        Parameters:
            pop (list): Current population, represented as a list
                of solutions (arrays of shape (L, 3)).
            scores (:obj:`np.ndarray`): Fitness functions associated
                to the individuals. Array thus has a length equal to
                the population size.

        Returns:
            :obj:`np.ndarray`: New mutated solution
                (array of shape (L, 3)).
        """
        # Shuffle the population
        indices = np.arange(len(pop))
        random.shuffle(indices)
        scores = scores[indices]

        # Elect a winner in each of the two partitions
        ps = self.partition_size
        left_winner = pop[indices[np.argmax(scores[:ps])]]
        right_winner = pop[indices[ps + np.argmax(scores[ps:2*ps])]]

        # Apply the cross-over and mutation operators
        individual = self.cross_over(left_winner, right_winner)
        return self.mutate(individual)

    def run(self, model, verbose=True):
        """Run heuristic optimizer on an initial solution,
        with given objective function.

        Parameters:
            model (:obj:`gaussfold.Model`): Gaussian model
            verbose (bool): Whether to display messages in stdout.

        Returns:
            :obj:`np.ndarray`: Optimal solution.
        """
        # Randomly initializes population and adds initial
        # solution to it
        initial_solution = model.get_coords()
        pop = [self.random_sol(initial_solution) for i in range(self.pop_size-1)]
        pop += [initial_solution]
        obj = model.evaluate

        # Set initial solution as the best one so far
        self.scores = list()
        best_score = -np.inf
        best_iteration = 0


        # Compute fitness functions on all individuals
        scores = np.asarray([obj(ind) for ind in pop])

        for k in range(self.n_iter):
            # Create new solution to replace worst solution
            new_ind = self.new_sol(pop, scores)
            worst = np.argmin(scores)
            pop[worst][:] = new_ind
            scores[worst] = obj(new_ind)
            assert(not np.isnan(scores[worst]))

            # Check if improvement
            if scores[worst] > best_score:
                best_score = scores[worst]
                best_iteration = k
            if verbose and (k + 1) % 100 == 0:
                print('Log-likelihood at iteration %i: %f' \
                    % (k + 1, best_score))
            self.scores.append(best_score)

            if np.isnan(best_score):
                if verbose:
                    print('[Warning] Invalid value encountered in heuristic solver')
                break

            # Stop algorithm if no more improvement
            if k - best_iteration >= self.early_stopping:
                break

        # Fine-tune solution with L-BFGS
        best_coords = pop[np.argmax(scores)]
        if self.use_lbfgs:
            new_coords = lbfgs(best_coords, model, verbose=verbose)
            if obj(new_coords) > best_score:
                best_coords = new_coords

        # Update coordinates in model
        model.set_coords(best_coords)
