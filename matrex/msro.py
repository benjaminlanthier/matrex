"""
A set of efficient, matrix-based tools for calculating the reordering of a
matrix's rows. As much as possible, matrix-based computations are used instead
of loops.

The algorithm is iterative. At each step, it chooses the next index that is
part of the new ordering from a list of candidates. For efficiency, the
algorithm computes the function to score the candidates in parallel. The only
loop is the iterative `while` loop which builds the new order one-by-one.


Date created: 2023-06-16
"""

# pylint: disable=C0103

from operator import itemgetter
import numpy as np
import networkx as nx

try:
    import matplotlib.pyplot as plt

    plot = True
except ImportError:
    plot = False


def pseudodiameter(graph: nx.Graph, seed: int = 111) -> tuple[int, int]:
    """
    Purpose
    -------
    Find nodes `i` and `j` which have a distance close to the diameter. Based
    on [Mathematica](https://reference.wolfram.com/language/GraphUtilities/ref/PseudoDiameter.html)'s
    approach.

    Arguments
    ---------
    graph : nx.Graph
        The row graph of the matrix.

    seed : int
        The seed for reproducibility of the initial random node choice.

    Returns
    -------
    (start, farthestNode) : tuple
        The node indices that are separated by the pseudodiameter.
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
            finalDegrees = dict(graph.degree(nbunch=distances.keys()))
            start = min(finalDegrees, key=finalDegrees.get)  # type: ignore
            final_loop = True
        else:
            start = farthest_node
            distance = new_distance

    return (start, farthest_node)


def find_active_rows(reordered: np.ndarray, row_graph: np.ndarray) -> np.ndarray:
    """
    Purpose
    -------
    Find active row indices by analyzing the submatrix of the row_graph
    which corresponds to the current front.

    If an entry `(i,j)` of the row_graph is nonzero, this means row `i` is
    connected to row `j` (they have at least one shared index that has a
    nonzero entry in both rows). The function finds all such nonzero.

    Arguments
    ---------
    reordered : list
        The current rows that are part of the front.

    row_graph : np.ndarray
        The row graph of the matrix.

    Returns
    -------
    active : np.ndarray of ints
        The active nodes indices.
    """
    # Use the row graph to find indices of active nodes that aren't already part of the front.
    candidates = np.nonzero(np.any(row_graph[reordered], axis=0))[0]
    active = np.setdiff1d(ar1=candidates, ar2=reordered)
    return active


def change_in_frontsize(
    front: np.ndarray, candidate_rows: np.ndarray, remaining_rows: np.ndarray
) -> np.ndarray:
    """
    Purpose
    -------
    Calculate the change in row and column frontsizes. The data
    structure encodes the remaining rows of the matrix for each
    choice of current row.

    Parameters
    ----------
    front : np.ndarray
        The front matrix before adding the current row.

    candidate_rows : np.ndarray of np.ndarrays
        The new rows to consider.

    remaining_rows : np.ndarray of np.ndarrays
        The remaining rows.

    Returns
    -------
    rcgain : np.ndarray
        The change in front size corresponding to each of the
        active rows if they were to be assembled next in the front.

    Example
    -------
    >>> matrix = np.array([[0, 1, 1, 0, 0],
    >>>                    [1, 0, 0, 1, 0],
    >>>                    [1, 0, 0, 1, 1],
    >>>                    [1, 1, 1, 1, 0],
    >>>                    [0, 1, 0, 1, 0],
    >>>                    [1, 0, 1, 1, 1]])

    >>> front = [[0, 1, 1, 0, 0],
    >>>          [1, 0, 0, 1, 0]]

    >>> candidate_rows = [[1, 0, 0, 1, 1],
    >>>                   [1, 1, 1, 1, 0],
    >>>                   [0, 1, 0, 1, 0],
    >>>                   [1, 0, 1, 1, 1]]

    >>> remaining_rows = [[[1, 1, 1, 1, 0],
    >>>                    [0, 1, 0, 1, 0],
    >>>                    [1, 0, 1, 1, 1]],
    >>>                   [[1, 0, 0, 1, 1],
    >>>                    [0, 1, 0, 1, 0],
    >>>                    [1, 0, 1, 1, 1]],
    >>>                   [[1, 0, 0, 1, 1],
    >>>                    [1, 1, 1, 1, 0],
    >>>                    [1, 0, 1, 1, 1]],
    >>>                   [[1, 0, 0, 1, 1],
    >>>                    [1, 1, 1, 1, 0],
    >>>                    [0, 1, 0, 1, 0]]]

    * Note that each array in remaining_rows contains the remaining
    rows in the matrix after we remove the front and the row that
    corresponds to the same row in candidate_rows.

    * We can then use matrix multiplication to figure out which
    columns are already part of the front and which ones are new,
    as well as the columns which are now fully summed.
    """
    front_columns = np.any(a=front, axis=0)
    new_front_columns = candidate_rows.astype(int) @ (~front_columns).astype(int).T
    zero_columns_remaining = (~np.any(a=remaining_rows, axis=1)).astype(int)
    fully_summed = np.sum(a=zero_columns_remaining * candidate_rows, axis=1)
    rcgain = 1 + new_front_columns - 2 * fully_summed
    return rcgain


def find_front_columns(front, candidate_rows):
    """
    Purpose
    -------
    For a given row, calculate nonzero indices that already describe the
    nonzero columns in the front.

    Arguments
    ---------
    front : np.ndarray
        The current matrix that describes the front.

    candidate_rows : np.ndarray of np.ndarrays
        A multidimensional array holding the candidate rows.

    Returns
    -------
    front_columns : np.ndarray of ints
        The number of columns already in the front.

    Note
    ----
    The Boolean matrix is converted to integers because matrix multiplication
    reflects the logical operation that should be done to get the number of
    columns that have a non-zero value in the `candidate_rows` that are already
    in the front.
    """
    current_front_columns = np.any(a=front, axis=0)
    front_columns = candidate_rows.astype(int) @ current_front_columns.astype(int).T
    return front_columns


def calculate_ordering(
    graph: nx.Graph,
    matrix: np.ndarray,
    weights: list = [2, 1, 0.2],
    seed: int = 111,
    verbose: bool = False,
) -> list:
    """
    Purpose
    -------
    For a connected graph, calculate an ordering of the nodes.

    Arguments
    ---------
    graph : nx.Graph
        The graph to reorder.

    matrix : np.ndarray
        The matrix for the entire original graph.

    weights : list[float, float, float]
        The weights (positive values) in the priority function of the algorithm.

    seed : integer
        For reproducibility of the pseudodiameter.

    verbose : bool
        If True, prints out information about the process.

    Returns
    -------
    order : list of ints
        A permutation of the graph nodes (the rows of the input matrix).
    """
    W1, W2, W3 = weights
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

        # This deletes particular columns from each array so that
        # we have only the unordered indices left
        all_unordered_indices = np.setdiff1d(ar1=all_nodes, ar2=order)

        # Find elements of each row to remove
        remove_mask = (
            (all_unordered_indices == active[:, None]).any(axis=0).nonzero()[0]
        )

        # We will remove one element per row of this array
        unordered_indices = np.tile(A=all_unordered_indices, reps=(len(active), 1))
        # Keep all elements that are NOT part of the removeMask and then reshape
        # From: https://stackoverflow.com/a/36502927
        m, n = unordered_indices.shape
        keep_mask = np.arange(n) != remove_mask[:, None]
        unordered_indices = unordered_indices[keep_mask].reshape(m, -1)

        # Build the rows and remaining matrices for each row of the orderedMask
        candidate_rows = matrix[active]
        remaining_rows = matrix[unordered_indices]

        # Calculate function values in matricized form
        values = find_front_columns(front=front, candidate_rows=candidate_rows)
        deltas = change_in_frontsize(
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


def msro(
    input_matrix: np.ndarray,
    weights: list = [2, 1, 0.2],
    show_row_graph: bool = False,
    seed: int = 111,
    verbose: bool = False,
) -> list:
    """
    Purpose
    -------
    Apply the MSRO (Modified Sloan Row Ordering) algorithm on the input matrix `m`
    in order to find the row ordering that minimizes the row front size.

    Arguments
    ---------
    input_matrix : np.ndarray
        The matrix to reorder.

    weights : list[float, float, float]
        The weights (positive values) in the priority function of the algorithm.

    show_row_graph : bool
        If True, shows the row graph generated by `networkx`.

    seed : integer
        For reproducibility of the pseudodiameter.

    verbose : bool
        If True, prints out information about the process.

    Returns
    -------
    total_reordering : list
        The new ordering of the rows.

    Example
    -------
    >>> import numpy as np
    >>> from matrex import msro
    >>> m = np.array([[1, 0, 1, 1, 0, 0],
    >>>               [0, 1, 0, 1, 1, 0],
    >>>               [1, 0, 1, 1, 0, 1],
    >>>               [0, 1, 0, 0, 0, 0],
    >>>               [0, 0, 0, 1, 1, 1],
    >>>               [0, 0, 0, 0, 0, 1]])
    >>> new_order = msro(m)
    >>> reordered_m = m[new_order]
    >>> print(f"The new rows ordering : {new_order}")
    >>> print(f"The reordered matrix : {m}")

    * The output should be :

    >>> The new rows ordering : [3, 1, 4, 5, 2, 0]
    >>> The reordered matrix :
    >>> array([[0 1 0 0 0 0],
    >>>        [0 1 0 1 1 0],
    >>>        [0 0 0 1 1 1],
    >>>        [0 0 0 0 0 1],
    >>>        [1 0 1 1 0 1],
    >>>        [1 0 1 1 0 0]])
    """
    # Only consider the nonzero pattern of elements
    matrix = input_matrix.astype(bool)

    # Compute initial distances between rows of the matrix.
    row_graph = nx.from_numpy_array(matrix @ matrix.T)
    if plot and show_row_graph:
        nx.draw(row_graph, with_labels=True)
        plt.show()

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
            graph=graph, matrix=matrix, weights=weights, seed=seed, verbose=verbose
        )
        total_reordering.extend(order)
    return total_reordering
