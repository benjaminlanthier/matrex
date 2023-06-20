"""
Purpose: Demonstrate an example of the reordering algorithm for a random matrix.

Date created: 2023-06-16
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import time
import networkx as nx
from utilities import MSRO
from vectorizedUtilities import vectorizedMSRO


def getWidths(matrix):
    widths = []
    for row in matrix:
        elements = np.nonzero(row)[0]
        if len(elements) > 0:
            width = np.max(elements) - np.min(elements)
        else:
            width = 0
        widths.append(width)
    return widths


if __name__ == "__main__":
    seed = 111  # For reproducibility of the pseudodiameter
    N = 50
    M = 50
    np.random.seed(seed=seed)
    matrix = np.zeros((M, N), dtype=int)
    for i in range(M):
        indices = np.random.choice(a=range(N), size=3, replace=False)
        matrix[i, indices] = 1

    """
    # Small test matrix
    matrix = np.array(
        [
            [1, 0, 1, 1, 0, 0],
            [0, 1, 0, 1, 1, 0],
            [1, 0, 1, 1, 0, 1],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 1, 1],
            [0, 0, 0, 0, 0, 1],
        ]
    )"""

    W1, W2, W3 = 2, 1, 0.2
    algorithmTime = time.time()
    orderSequential = MSRO(input_matrix=matrix.T, W1=W1, W2=W2, W3=W3, seed=seed)
    order = vectorizedMSRO(input_matrix=matrix.T, W1=W1, W2=W2, W3=W3, seed=seed)
    print("Order: ", order)
    print("Seque: ", orderSequential)
    print(
        "Sequential and vectorized orderings the same? ",
        np.allclose(order, orderSequential),
    )
    finishTime = time.time()
    print("Total algorithm time: {}s".format(finishTime - algorithmTime))

    # Create a figure and post results
    # Original stats
    widths = getWidths(matrix)

    # Reduced fill-in stats
    reordered = matrix.T[order].T
    reorderedWidths = getWidths(reordered)

    print("Original matrix average width : ", np.mean(widths))
    print("Reordered matrix average width: ", np.mean(reorderedWidths))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8))
    ax1.imshow(matrix, cmap="binary")
    ax2.imshow(reordered, cmap="binary")
    ax1.set_title(r"Original matrix ($W = {:.2f}$)".format(np.mean(widths)))
    ax2.set_title(r"Reduced matrix ($W = {:.2f}$)".format(np.mean(reorderedWidths)))
    plt.savefig("reorderedMatrixN{}.pdf".format(N), dpi=300)
    # plt.show()
