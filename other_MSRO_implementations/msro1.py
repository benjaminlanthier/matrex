"""
Purpose: A set of efficient, matrix-based tools for calculating the reordering of a
         matrix's rows. As much as possible, I use matrix-based computations instead
         of loops.
Date created: 2023-06-14
"""

# pylint: disable=C0103

import numpy as np
import networkx as nx


def change_in_frontsize(front, row, remaining):
    """
    - Purpose: Calculate the change in row and column frontsizes.
    - Input:
      - front (array): The front matrix before adding the current row.
      - row (array): The new row to consider.
      - remaining (array): The remaining rows.
    - Output:
      - (integer): The change in front sizes.
    """
    # Boolean mask for which indices are present
    front_columns = np.any(a=front, axis=0)

    new_columns = np.sum(np.logical_and(~front_columns, row))
    fully_summed = np.sum(np.logical_and(~np.any(a=remaining, axis=0), row))
    return 1 + new_columns - 2 * fully_summed


def pseudodiameter(graph, seed=111):
    """
    - Purpose: Find nodes i and j which have a distance close to the diameter.
               Based on Mathematica's approach:
               https://reference.wolfram.com/language/GraphUtilities/ref/PseudoDiameter.html
    - Input:
      - graph (NetworkX object): The row graph of the matrix.
      - seed (integer): The seed for reproducibility of the initial random node choice.
    - Output:
      - (i, j) (tuple): The node indices that are separated by the pseudodiameter.
    """
    np.random.seed(seed=seed)
    nodes = graph.nodes()
    start = np.random.choice(nodes)
    final_loop = False
    distance = 0
    while True:
        distances = nx.algorithms.single_source_shortest_path_length(
            G=graph, source=start
        )
        farthest_node = max(distances, key=distances.get)
        new_distance = distances[farthest_node]
        if final_loop:
            break
        if new_distance == distance:
            final_degrees = dict(graph.degree(nbunch=distances.keys()))
            start = min(final_degrees, key=final_degrees.get)
            final_loop = True
        else:
            start = farthest_node
            distance = new_distance

    # print(f"(start, target) : ({start}, {farthestNode})")

    return (start, farthest_node)


def find_active_rows(reordered, row_graph):
    """
    - Purpose: Find active row indices by analyzing the submatrix of the row_graph
               which corresponds to the current front.

               If an entry (i,j) of the row_graph is nonzero, this means row i is
               connected to row j (they have at least one shared index that has a
               nonzero entry in both rows). First, the function finds all such nonzero
    - Input:
      - reordered (list of integers): The current nodes that are part of the front.
      - row_graph (array): The row graph of the matrix.
    - Output:
      - active (list of integers): The active node indices.
    """
    # Use the row graph to find indices of active nodes that aren't already part of the front.

    candidates = np.nonzero(np.any(row_graph[reordered], axis=0))[0]
    active = np.setdiff1d(ar1=candidates, ar2=reordered)
    return active


def find_front_columns(matrix, row):
    """
    - Purpose: For a given row, calculate the number of columns that
               have nonzero elements and are already in the front.
    - Input:
      - matrix (array): The current matrix that describes the front.
      - row (array): The particular row to consider.
    - Output:
      - (integer): The number of columns already in the front.
    """
    current_front_columns = np.any(a=matrix, axis=0)  # Boolean mask
    return np.sum(np.logical_and(current_front_columns, row))


def calculate_ordering(graph, matrix, W1=2, W2=1, W3=0.2, seed=111, verbose=False):
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
        priorities = {}
        front = matrix[order]
        for node in active:
            order.append(node)
            unordered = np.setdiff1d(ar1=all_nodes, ar2=order)
            order.pop()
            row, remaining = matrix[node], matrix[unordered]
            front_columns = find_front_columns(matrix=front, row=row)
            delta_front = change_in_frontsize(front=front, row=row, remaining=remaining)
            score = -W1 * delta_front + W2 * distances[node] - W3 * front_columns
            priorities[node] = score
            if verbose:
                print("Node: ", node)
                print("Distance to target: ", distances[node])
                print(
                    "Columns already in the front: ",
                    find_front_columns(matrix=front, row=row),
                )
                print("Order    : ", order)
                print("Unordered: ", unordered)
        # print(f"In sequential : {np.array(data)}")
        # Choose highest priority
        selection = max(priorities, key=priorities.get)
        order.append(selection)
        if verbose:
            print("Start, target: ", start, target)
            print("Active: ", active)
            print("Priorities: ")
            print(priorities)
            print("Selection: ", selection)
            print()
    assert len(np.unique(order)) == len(all_nodes)
    return order


def msro(input_matrix, W1=2, W2=1, W3=0.2, seed=111, verbose=False):
    """
    - Purpose: Reorder the rows of the input matrix to minimize its front.
    - Input:
      - input_matrix (array): The matrix to reorder.
      - W1, W2, W3 (positive integers): The weights for balancing the algorithm's
                                        considerations.
      - seed (integer): For reproducibility of the pseudodiameter.
      - verbose (Boolean): If True, prints out information about the process.
    - Output:
      - totalReordering (1D array): The permutation of the rows.
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
        order = calculate_ordering(
            graph=graph, matrix=matrix, W1=W1, W2=W2, W3=W3, seed=seed, verbose=verbose
        )
        total_reordering.extend(order)
    return total_reordering
