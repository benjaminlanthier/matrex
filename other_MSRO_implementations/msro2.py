"""
A set of efficient, matrix-based tools for calculating the reordering of a
matrix's rows. As much as possible, matrix-based computations are used instead
of loops.

The algorithm is iterative. At each step, it chooses the next index that is
part of the new ordering from a list of candidates. This implementation follows
an object that is a list of lists where the lists contain 2 values. The first
one is the original row index and the second one is its state, either 'inactive',
'active' or 'assembled'. It moves those sublists in the list so that the index of
the sublist corresponds to the physical (current) row index in the matrix. It lets
the algorithm follow the original row index wherever this row is in the matrix.
On my laptop, it deals with a 1000x1000 matrix in around 4 minutes.


Date created: 2023-06-16
"""

# pylint: disable=C0103

from operator import itemgetter
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from msro1 import pseudodiameter


def generate_rowGraph(m: np.ndarray, show: bool = False) -> nx.Graph:
    """
    Purpose
    -------
    Generate the row grpah and extract each of its subgraphs in a list.

    Parameters
    ----------
    m : np.ndarray
        The input matrix on which the MSRO algorithm will be applied.

    show : bool
        To either show the row graph (`True`) of not (`False`).

    Returns
    -------
    subgraphs : list of nx.Graph objects
        The list of all the subgraphs of the row graph.
    """
    edges_list = np.argwhere(np.triu(m @ m.T, k=1))
    Gr = nx.Graph()
    Gr.add_nodes_from(range(m.shape[0]))
    Gr.add_edges_from(edges_list)
    if show:
        nx.draw(Gr, with_labels=True)
        plt.show()
    sub_graphs = [
        Gr.subgraph(component).copy() for component in nx.connected_components(Gr)
    ]
    return sub_graphs


def get_subgraph_data(
    subgraph: nx.Graph, seed: int = 111
) -> tuple[int, int, dict, dict]:
    """
    Purpose
    -------
    Extract important information about the row graph of the input matrix
    so that it can be used in the algorithm.

    Parameters
    ----------
    subgraph : nx.Graph
        The subgraph of the entire row graph of `m` on which the MSRO algorithm
        will be applied. Basically, if the row graph is disconnected, it applied the
        algorithm on each subgraphs of the row graph.

    silent : bool
        Either show the time needed for each steps in the function (True)
        or not (False).

    Returns
    -------
    s, e : int, int
        The start (s) and target (e) nodes selected for the MSRO algorithm.

    distances : dict
        The distances between each nodes of the subgraph and the target node e.

    neighbors : dict
        The neighbors of all the nodes in the subgraph of the row graph.
    """
    s, e = pseudodiameter(subgraph, seed=seed)
    # Generate the dictionnary of the distances of every nodes from the target node
    distances = nx.shortest_path_length(subgraph, target=e)
    neighbors = {}
    for node in subgraph.nodes():
        neighbors_i = nx.neighbors(subgraph, node)
        neighbors[node] = list(neighbors_i)

    return s, e, distances, neighbors


def update_rows_order(
    rows_states: np.ndarray,
    row_to_assemble_id: int,  # the current row id, not the original one
    assembling_step: int,
    neighbors: dict,
) -> np.ndarray:
    """
    Purpose
    -------
    Update the `rows_states` object so that the original row indexes (and
    their states) are saved with the right current row index.

    Parameters
    ----------
    rows_states : np.ndarray([int, str], dtype = object)
        The object that saves the original row indexes and their states.

    row_to_assemble_id : int
        The next row to assemble in the matrix on which we want to apply the
        MSRO algorithm.

    assembling_step : int
        The assembling step the algorithm is at.

    neighbors : dict
        The neighbors of each nodes in the row graph.

    Returns
    -------
    rows_states_copy : np.ndarray([int, str], dtype = object)
        The updated `rows_states` object where we assembled the row with index
        `row_to_assembled_id` and updated its neighbors to have the state "active".
    """
    # Update the original rows ids in the current rows ids
    original_row_id = rows_states[row_to_assemble_id, 0]
    if rows_states[assembling_step, 0] != original_row_id:
        rows_states[assembling_step + 1 : row_to_assemble_id + 1] = rows_states[
            assembling_step:row_to_assemble_id
        ]
        rows_states[assembling_step, 0] = original_row_id

    # Update the original rows states (inactive, active or assembled)
    rows_states[assembling_step][1] = "assembled"
    current_rows_ids = np.where(
        np.isin(rows_states[:, 0], neighbors[original_row_id]) == True
    )[0]
    not_assembled_mask = rows_states[current_rows_ids, 1] != "assembled"
    rows_states[current_rows_ids[not_assembled_mask], 1] = "active"

    return rows_states


def get_rcgain_and_nold(
    m: np.ndarray,
    rows_states: np.ndarray,
    row_ids: list,  # the current row indexes
    nb_assembled_rows: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Purpose
    -------
    Evaluate the increase of the front row size and the column row size from the
    permutation of the 2 given rows of the matrix.

    Parameters
    ----------
    `m` : np.ndarray
        The given matrix for which we want to evaluate the increase in the front
        row and column sizes if we were to permute the 2 given rows.

    rows_states : np.ndarray([int, str], dtype = object)
        The object that saves the original row indexes and their states.

    nb_assembled_rows : int
        The number of rows that have already been assembled in the front.

    Returns
    -------
    `rcgain` : list
        The increase of the front row size and the column row size from the
        assembling of all the next active rows.

    `nold` : list
        The number of variables in row `i` that are candidates for elimination
        and that are already in the front.
    """
    new_rows_order = [i[0] for i in rows_states]
    new_orders = []
    tensor = []
    for row_id in row_ids:
        new_order = np.delete(new_rows_order, row_id)
        new_order = np.insert(new_order, nb_assembled_rows, rows_states[row_id, 0])
        new_orders.append(new_order)

    tensor = np.tile(A=m, reps=(len(row_ids), 1, 1))
    tensor = tensor[np.arange(tensor.shape[0])[:, None], new_orders]

    in_front = tensor[:, : nb_assembled_rows + 1]
    in_front_counts = np.sum(in_front, axis=1)
    newc = np.count_nonzero(
        (tensor[:, nb_assembled_rows] == 1) & (in_front_counts == 1), axis=1
    )

    not_in_front = tensor[:, nb_assembled_rows + 1 :]
    not_in_front_counts = np.sum(not_in_front, axis=1)
    s = np.count_nonzero(
        (tensor[:, nb_assembled_rows] == 1) & (not_in_front_counts == 0), axis=1
    )

    rgain = 1 - s
    cgain = newc - s
    rcgain = rgain + cgain

    new_rows = tensor[:, nb_assembled_rows]
    has_ones = np.logical_and(new_rows[:, None], tensor[:, :nb_assembled_rows]).any(
        axis=1
    )
    nold = np.sum(has_ones, axis=1)

    return rcgain, nold


def initialize_P(
    m: np.ndarray, rows_states: np.ndarray, distances: dict, weights: list, count: int
) -> tuple[list, float]:
    """
    Purpose
    -------
    Evaluate the initial values of P (priority function) for all the rows of `m`.

    Parameters
    ----------
    m : np.ndarray
        The matrix on which we want to apply the MSRO algorithm.

    rows_states : np.ndarray([int, str], dtype = object)
        The object that saves the original row indexes and their states.

    neighbors : dict
        The neighbors of each nodes in the row graph.

    distances : dict
        The distances between every nodes in the row graph and the target node `e`.

    weights : list
        The weights `W1`, `W2` and `W3` in the priority function.

    count : int
        Value to keep track of the assembling step through the subgraphs of the
        row graph.

    Returns
    -------
    P : np.ndarray
        The initial priority function values.

    v : float
        The normalizing factor in the expression of the priority function.
    """
    W1, W2, W3 = weights
    P = []
    v = 1
    P = np.array([np.NINF] * m.shape[0])
    current_rows_ids = np.where(np.isin(rows_states[:, 0], list(distances.keys())))[0]

    original_rows_ids = rows_states[current_rows_ids, 0]
    rcgain, nold = get_rcgain_and_nold(m, rows_states, current_rows_ids, count)
    P[original_rows_ids.astype(int)] = (
        -W1 * rcgain + W2 * v * itemgetter(*original_rows_ids)(distances) - W3 * nold
    )
    return P, v


def get_row_to_assemble(
    P: list, rows_states: np.ndarray, s: int, assembling_step: int, count: int
) -> int:
    """
    Purpose
    -------
    Compute the next row to assemble in the front according to the priority
    function.

    Parameters
    ----------
    P : list
        The values of the priority function for each rows in the matrix.

    rows_states : np.ndarray([int, str], dtype = object)
        The object that saves the original row indexes and their states.

    s : int
        The start node in the row graph.

    assembling_step : int
        The assembling step the algorithm is at.

    count : int
        Value to keep track of the assembling step through the subgraphs of the
        row graph.

    Returns
    -------
    current_row_to_assemble : int
        The next row index to assemble in the front.
    """
    if assembling_step == count:
        current_row_to_assemble = np.where(rows_states[:, 0] == s)[0][0]
    else:
        indices = np.where(rows_states[:, 1] != "active")[0]
        not_active_rows_idx = rows_states[indices, 0]
        P_copy = np.array(P, copy=True)
        P_copy[not_active_rows_idx.astype(int)] = np.NINF
        original_row_to_assemble = np.argmax(P_copy)
        current_row_to_assemble = np.where(
            rows_states[:, 0] == original_row_to_assemble
        )[0][0]
    return current_row_to_assemble


def update_P(
    m: np.ndarray,
    P: list,
    v: float,
    rows_states: np.ndarray,
    weights: list,
    distances: dict,
) -> list:
    """
    Purpose
    -------
    Update the values of P (priority function) for all the active rows in `m`.
    Specify the value `np.NINF` to an assembled row.

    Parameters
    ----------
    m : np.ndarray
        The matrix on which we want to apply the MSRO algorithm.

    P : list
        The values of the priority function for each rows in the matrix.

    v : float
        The normalizing factor in the expression of the priority function.

    rows_states : np.ndarray([int, str], dtype = object)
        The object that saves the original row indexes and their states.

    neighbors : dict
        The neighbors of each nodes in the row graph.

    weights : list
        The weights `W1`, `W2` and `W3` in the priority function.

    distances : dict
        The distances between every nodes in the row graph and the target node `e`.

    Returns
    -------
    P : list
        The updated priority function values.
    """
    indices = np.where(rows_states[:, 1] == "assembled")[0]
    original_assembled_rows_idx = rows_states[indices, 0]
    nb_assembled_rows = len(original_assembled_rows_idx)
    last_assembled_row_idx = original_assembled_rows_idx[-1]
    P[last_assembled_row_idx] = np.NINF
    current_active_rows_idx = np.where(rows_states[:, 1] == "active")[0]
    W1, W2, W3 = weights
    original_active_rows_idx = rows_states[current_active_rows_idx, 0]
    rcgain, nold = get_rcgain_and_nold(
        m, rows_states, current_active_rows_idx, nb_assembled_rows
    )
    P[original_active_rows_idx.astype(int)] = (
        -W1 * rcgain
        + W2 * v * itemgetter(*original_active_rows_idx)(distances)
        - W3 * nold
    )
    return P


def msro(
    m: np.ndarray,
    perm: str = "rows",
    weights: list = [2, 1, 0.2],
    seed: int = 111,
    show_Gr: bool = False,
    verbose: bool = False,
) -> list:
    """
    Purpose
    -------
    Apply the MSRO (Modified Sloan Row Ordering) algorithm on the input matrix `m`
    in order to find the row (or column) ordering that minimizes the row (or column)
    front size.

    Parameters
    ----------
    m : np.ndarray
        The matrix on which we want to apply the MSRO algorithm.

    perm : str
        Can be either `"rows"` or `"columns`". It lets the user specify if
        the rows or the columns of `m` should be reordered. If `"columns"` is
        specified, then the algorithm is applied on `m.T`.

    weights : list
        The weights `W1`, `W2` and `W3` in the priority function. Values are set
        by default to be 2, 1 and 0.2 for `W1`, `W2` and `W3` respectively.

    seed : int
        The seed for the random node selection in `pseudodiameter()`.

    show_Gr : bool
        The either show the row graph (True) or not (False).

    verbose : bool
        Either to print information on the terminal (True) or not (False).

    Returns
    -------
    final_order : list
        The list of reordered rows (or columns) for `m`.

    Example
    -------
    >>> import numpy as np
    >>> from matrexalgs import msro
    >>> m = np.array([[1, 0, 1, 1, 0, 0],
    >>>               [0, 1, 0, 1, 1, 0],
    >>>               [1, 0, 1, 1, 0, 1],
    >>>               [0, 1, 0, 0, 0, 0],
    >>>               [0, 0, 0, 1, 1, 1],
    >>>               [0, 0, 0, 0, 0, 1]])
    >>> new_order = msro(m, perm = "rows")
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
    if perm == "columns":
        matrix = np.array(m.T, dtype=bool, copy=True)
    elif perm == "rows":
        matrix = np.array(m, dtype=bool, copy=True)
    else:
        raise f"Bad input for `perm`. It should be either 'rows' or 'columns', but you used {perm}"

    rows_states = []
    for i in range(matrix.shape[0]):
        rows_states.append([i, "inactive"])
    rows_states = np.array(rows_states, dtype=object)

    Gr_subgraphs = generate_rowGraph(matrix, show=show_Gr)
    count = 0
    for subgraph in Gr_subgraphs:
        s, e, distances, neighbors = get_subgraph_data(subgraph, seed=seed)
        for assembling_step in range(count, count + len(subgraph)):
            if assembling_step == count:
                P, v = initialize_P(matrix, rows_states, distances, weights, count)
                row_id = get_row_to_assemble(P, rows_states, s, assembling_step, count)
                rows_states = update_rows_order(
                    rows_states, row_id, assembling_step, neighbors
                )
                if verbose:
                    print(f"-------Step {assembling_step}-------")
                    print(f"Initial P : {P}")
                    print(f"Next row to assemble : {row_id}")
                    print(f"Next rows order : \n{rows_states}")
            else:
                P = update_P(matrix, P, v, rows_states, weights, distances)
                row_id = get_row_to_assemble(P, rows_states, s, assembling_step, count)
                rows_states = update_rows_order(
                    rows_states, row_id, assembling_step, neighbors
                )
                if verbose:
                    print(f"\n-------Step {assembling_step}-------")
                    print(f"Updated P : {P}")
                    print(f"Next row to assemble : {row_id}")
                    print(f"Next rows order : \n{rows_states}")
        count += len(subgraph)

    final_order = [v[0] for v in rows_states]
    return final_order


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
