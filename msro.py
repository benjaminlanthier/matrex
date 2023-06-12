"""
Module that groups matrix reordering algorithms. Right now, it contains the [MSRO](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1506(199904/05)6:3%3C189::AID-NLA160%3E3.0.CO;2-C)
algorithm, which was first implemented in the [HSL](https://www.hsl.rl.ac.uk/catalogue/mc62.html)
library.

`MatRexAlg` is for Matrix Reordering Algorithms.

By - Benjamin Lanthier
"""

# pylint: disable=C0103

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from msro_utils import (initialize_rowGraph,
                        get_nodes_s_and_e,
                        get_row_to_assemble,
                        initialize_P,
                        update_data)


def get_mean_row_front_size(m):
    """
    Purpose
    -------
    Get the mean row front size of the given matrix.

    Parameters
    ----------
    `m` : np.ndarray
        The input matrix for which we want to evaluate the mean front row size.

    Returns
    -------
    np.mean(front_row_sizes) : float
        The mean of the front row size of `m`.
    """
    front_row_sizes = []
    for row in m:
        ones_indices = np.where(row == 1)[0]
        front_row_sizes.append(ones_indices[-1] - ones_indices[0])
    return np.mean(front_row_sizes)


def get_max_row_front_size(m):
    """
    Purpose
    -------
    Get the max row front size of the given matrix.

    Parameters
    ----------
    `m` : np.ndarray
        The input matrix for which we want to evaluate the mean front row size.

    Returns
    -------
    np.max(front_row_sizes) : float
        The max of the front row size of `m`.
    """
    front_row_sizes = []
    for row in m:
        ones_indices = np.where(row == 1)[0]
        front_row_sizes.append(ones_indices[-1] - ones_indices[0])
    return np.max(front_row_sizes)


def msro(
    m: np.ndarray,
    perm: str = "columns",
    show_rowGraph: bool = False,
    smallest_value: int = -1e12,
) -> np.ndarray:
    """
    Purpose
    -------
    Apply the [MSRO](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1506(199904/05)6:3%3C189::AID-NLA160%3E3.0.CO;2-C)
    algorithm on the given matrix `m` in order to minimize the row front size and
    the front column size at the same time. This algorithm can permute the rows or
    the columns of `m`.

    Parameters
    ----------
    `m` : np.ndarray
        The matrix on which we want to apply the MSRO algorithm.

    `perm` : str
        Whether to permute the rows or the columns of the given matrix. Set to
        `columns` by default for the sites relabeling in an MPS-MPO problem.

    `show_rowGraph` : bool
        Whether to show (True) the row graph or not (False).

    `smallest_value` : int
        The smallest value given to assembled nodes in the priodity function.

    Returns
    -------
    `m` : np.ndarray
        The modified version of `m` after the application of the algorithm.
    """
    # Properly setup the given matrix to permute the rows or the columns
    if perm == "columns":
        matrix = np.array(m, copy=True).T
    elif perm == "rows":
        matrix = np.array(m, copy=True)
    else:
        error_message = f"""'{perm}' is not an acceptable input for msro().\
                            Please choose either `rows` or `columns`"""
        raise ValueError(error_message)

    # Initialize the data for the algorithm
    Gr = initialize_rowGraph(matrix)
    if show_rowGraph:
        nx.draw(Gr, with_labels=True)
        plt.show()
    s, e = get_nodes_s_and_e(Gr)
    P = initialize_P(matrix, Gr, e)

    # Update those variables until P contains only -np.inf
    while not np.all(P == smallest_value):
        assembling_idx, row_to_assemble_idx = get_row_to_assemble(Gr, P, s)
        rows_to_switch = [assembling_idx, row_to_assemble_idx]
        matrix, P, Gr, e = update_data(
            matrix, P, Gr, rows_to_switch, e, smallest_value=smallest_value
        )

    # Return the matrix in the same format it was taken as an input
    if perm == "columns":
        return matrix.T
    elif perm == "rows":
        return matrix
