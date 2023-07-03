"""
Purpose: Demonstrate an example of the reordering algorithm for a random
matrix. It shows the form of the output and that all implementations give
the same results.

Date created: 2023-06-16

"""

# pylint: disable=C0103

import sys, os

sys.path.append(os.path.realpath(os.path.dirname(__file__) + "/.."))

import numpy as np
from msro1 import msro as msro1
from msro2 import msro as msro2

from matrex.msro import msro  # type: ignore


# # Kept it so it is still accessible
# def get_widths(matrix):
#     widths = []
#     for row in matrix:
#         elements = np.nonzero(row)[0]
#         if len(elements) > 0:
#             width = np.max(elements) - np.min(elements)
#         else:
#             width = 0
#         widths.append(width)
#     return widths


if __name__ == "__main__":
    seed = 666  # For reproducibility of the pseudodiameter
    np.random.seed(seed=seed)
    nb_simulations = 100
    min_size = 20
    max_size = 20  # if max_size == min_size, only one matrix size will be generated per simulation
    incrementation = 4

    weights = [2, 1, 0.2]
    W1, W2, W3 = weights
    for size in range(min_size, max_size + 1, incrementation):
        comparison = []
        for count in range(nb_simulations):
            matrix = np.zeros((size, size), dtype=int)
            for i in range(size):
                indices = np.random.choice(a=range(size), size=3, replace=False)
                matrix[i, indices] = 1
            # Theoretical toy example
            # matrix = np.array([[1, 0, 1, 1, 0, 0],
            #                    [0, 1, 0, 1, 1, 0],
            #                    [1, 0, 1, 1, 0, 1],
            #                    [0, 1, 0, 0, 0, 0],
            #                    [0, 0, 0, 1, 1, 1],
            #                    [0, 0, 0, 0, 0, 1]])
            msro1_order = msro1(input_matrix=matrix, W1=W1, W2=W2, W3=W3, seed=seed)
            msro2_order = msro2(matrix, perm="rows", seed=seed)
            msro_order = msro(input_matrix=matrix, weights=weights, seed=seed)
            comparison.append([msro_order == msro1_order, msro_order == msro2_order])
            if count % 100 == 0:
                print(f"msro1 order : {msro1_order}")
                print(f"msro2 order : {msro2_order}")
                print(f"msro order  : {msro_order}")
                print()

        results = np.sum(comparison, axis=0)
        print(f"Number of permutations that were the same (size = {size}) : {results}")
        print(f"100% of the permutations are the same : {results == nb_simulations}")
