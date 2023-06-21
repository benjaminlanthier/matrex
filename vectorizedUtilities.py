"""
Purpose: A set of efficient, matrix-based tools for calculating the reordering of a
         matrix's rows. As much as possible, I use matrix-based computations instead
         of loops.

         The algorithm is iterative. At each step, it chooses the next index that is
         part of the new ordering from a list of candidates. For efficiency, the
         algorithm computes the function to score the candidates in parallel. The only
         loop is the iterative while loop which builds the new order one-by-one.


Date created: 2023-06-16
"""

import numpy as np
import networkx as nx
from utilities import pseudodiameter, find_active_rows
from operator import itemgetter


def change_in_frontsize_vectorized(front, candidate_rows, remaining_rows):
    """
    - Purpose: Calculate the change in row and column frontsizes. The
               data structure encodes the remaining rows of the matrix
               for each choice of current row.

               For example:

               matrix =         [[0, 1, 1, 0, 0],
                                 [1, 0, 0, 1, 0],
                                 [1, 0, 0, 1, 1],
                                 [1, 1, 1, 1, 0],
                                 [0, 1, 0, 1, 0],
                                 [1, 0, 1, 1, 1]]

               front =          [[0, 1, 1, 0, 0],
                                 [1, 0, 0, 1, 0]]

               candidate_rows = [[1, 0, 0, 1, 1],
                                 [1, 1, 1, 1, 0],
                                 [0, 1, 0, 1, 0],
                                 [1, 0, 1, 1, 1]]

               remaining_rows = [
                                 [[1, 1, 1, 1, 0],
                                  [0, 1, 0, 1, 0],
                                  [1, 0, 1, 1, 1]],

                                 [[1, 0, 0, 1, 1],
                                  [0, 1, 0, 1, 0],
                                  [1, 0, 1, 1, 1]],

                                 [[1, 0, 0, 1, 1],
                                  [1, 1, 1, 1, 0],
                                  [1, 0, 1, 1, 1]],

                                 [[1, 0, 0, 1, 1],
                                  [1, 1, 1, 1, 0],
                                  [0, 1, 0, 1, 0]]
                                ]

                Note that each array in remaining_rows contains the remaining
                rows in the matrix after we remove the front and the row that
                corresponds to the same row in candidate_rows.

                We can then use matrix multiplication to figure out which
                columns are already part of the front and which ones are new,
                as well as the columns which are now fully summed.
    - Input:
      - front (array): The front matrix before adding the current row.
      - candidate_rows (array of arrays): The new rows to consider.
      - remaining_rows (array of arrays): The remaining rows.
    - Output:
      - (integer): The change in front size.
    """
    front_columns = np.any(a=front, axis=0)
    new_front_columns = candidate_rows.astype(int) @ (~front_columns).astype(int).T
    zero_columns_remaining = (~np.any(a=remaining_rows, axis=1)).astype(int)
    # print("Front columns shape: ", front_columns.shape)
    # print("New front columns  : ", new_front_columns.shape)
    # print("Zero columns shape : ", zero_columns_remaining.shape)
    # print("Multiplication: ")
    # print(zero_columns_remaining @ candidate_rows.astype(int).T)
    # fully_summed = np.sum(
    #    a=zero_columns_remaining @ candidate_rows.astype(int).T, axis=1
    # )
    # fully_summed = np.diag(
    #    v=zero_columns_remaining @ candidate_rows.astype(int).T
    # )
    # Is this the fastest?
    fully_summed = np.sum(a=zero_columns_remaining * candidate_rows, axis=1)
    # print("Fully summed: ")
    # print(fully_summed)

    return 1 + new_front_columns - 2 * fully_summed


def find_front_columns_vectorized(front, candidate_rows):
    """
    - Purpose: For a given row, calculate nonzero indices that already
               describe the nonzero columns in the front.

               Note: I convert the Boolean matrices to integers because
                     matrix multiplication reflects the logical operation
                     I want to do.
    - Input:
      - front (array): The current matrix that describes the front.
      - candidate_rows (array of arrays): A multidimensional array
                                          holding the candidate rows.
    - Output:
      - front_columns (array of integers): The number of columns already in the front.
    """
    current_front_columns = np.any(a=front, axis=0)
    return candidate_rows.astype(int) @ current_front_columns.astype(int).T


def calculate_ordering_vectorized(
    graph, matrix, W1=2, W2=1, W3=0.2, seed=111, verbose=False
):
    """
    - Purpose: For a connected graph, calculate an ordering of the nodes.
    - Inputs:
      - graph (NetworkX object): The graph to reorder.
      - matrix (array): The matrix for the entire original graph.
      - W1, W2, W3 (positive integers): The weights for balancing the algorithm's
                                        considerations.
      - seed (integer): For reproducibility of the pseudodiameter.
      - verbose (Boolean): If True, prints out information about the process.
    - Output:
      - order (list of integers): A permutation of the graph nodes.
    """
    row_graph = matrix @ matrix.T
    start, target = pseudodiameter(graph=graph, seed=seed)
    distances = nx.algorithms.single_source_shortest_path_length(G=graph, source=target)
    if verbose:
        print("Distances: ", distances)
    all_nodes = graph.nodes()
    order = [start]
    while len(order) < len(graph):
        # Find active rows
        active = find_active_rows(reordered=order, row_graph=row_graph)
        front = matrix[order]

        # Still need to understand these sorting lines!
        # This deletes particular columns from each array so that
        # we have only the unordered indices left
        all_unordered_indices = np.setdiff1d(ar1=all_nodes, ar2=order)

        # Find elements of each row to remove
        removeMask = (all_unordered_indices == active[:, None]).any(axis=0).nonzero()[0]

        # We will remove one element per row of this array
        unordered_indices = np.tile(A=all_unordered_indices, reps=(len(active), 1))
        # Keep all elements that are NOT part of the removeMask and then reshape
        # From: https://stackoverflow.com/a/36502927
        m, n = unordered_indices.shape
        keep_mask = np.arange(n) != removeMask[:, None]
        unordered_indices = unordered_indices[keep_mask].reshape(m, -1)

        # Build the rows and remaining matrices for each row of the orderedMask
        candidate_rows = matrix[active]
        remaining_rows = matrix[unordered_indices]

        # Calculate function values in matricized form
        values = find_front_columns_vectorized(
            front=front, candidate_rows=candidate_rows
        )
        deltas = change_in_frontsize_vectorized(
            front=front, candidate_rows=candidate_rows, remaining_rows=remaining_rows
        )
        score = -W1 * deltas + W2 * itemgetter(*tuple(active))(distances) - W3 * values
        selectionIndex = np.argmax(a=score)
        selection = active[selectionIndex]
        order.append(selection)
        if verbose:
            print("Start, target: ", start, target)
            print("Active: ", active)
            print("Selection: ", selection)
            print()
    assert len(np.unique(order)) == len(all_nodes)
    return order


def vectorizedMSRO(input_matrix, W1=2, W2=1, W3=0.2, seed=111, verbose=False):
    """
    - Purpose: Reorder the rows of the input matrix to minimize its front.
    - Input:
      - input_matrix (array): The matrix to reorder.
      - W1, W2, W3 (positive integers): The weights for balancing the algorithm's
                                        considerations.
      - seed (integer): For reproducibility of the pseudodiameter.
      - verbose (Boolean): If True, prints out information about the process.
    - Output:
      - total_reordering (1D array): The permutation of the rows.
    """

    # Only consider the nonzero pattern of elements
    matrix = input_matrix.astype(bool)

    # Compute initial distances between rows of the matrix.
    row_graph = nx.from_numpy_array(matrix @ matrix.T)

    # Check if the graph is disconnected and break into components
    total_reordering = []
    subgraphs = [
        row_graph.subgraph(component).copy()
        for component in nx.connected_components(row_graph)
    ]
    if verbose:
        print("Component sizes: ", [len(i) for i in subgraphs])
    for e, graph in enumerate(subgraphs):
        if verbose:
            print("Subgraph: ", e)
        order = calculate_ordering_vectorized(
            graph=graph, matrix=matrix, W1=W1, W2=W2, W3=W3, seed=seed, verbose=verbose
        )
        total_reordering.extend(order)
    return total_reordering
