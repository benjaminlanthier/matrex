"""
Module that groups matrix reordering algorithms. Right now, it contains the [MSRO](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1506(199904/05)6:3%3C189::AID-NLA160%3E3.0.CO;2-C)
algorithm, which was first implemented in the [HSL](https://www.hsl.rl.ac.uk/catalogue/mc62.html)
library.
"""

# pylint: disable=C0103

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def initialize_rowGraph(m: np.ndarray) -> nx.Graph:
    """
    Purpose
    -------
    Generate the row graph of the given matrix `m`. A row graph is a graph that uses
    each rows of the matrix as a vertex. The node `i` has an edge between itself and
    the node `j` if the non-zero values of the node `i` overlap the non-zero values
    of the node `j` in any of their columns.

    * Here, we use the fact that the row graph of `matrix` corresponds to the
    undirected graph of `B` if `B = M @ M.T` to optimize the generation of `Gr`.

    Parameters
    ----------
    `m` : np.ndarray
        The matrix for which we want to get the row-graph.

    Returns
    -------
    `Gr` : nx.Graph
        The row-graph of the given matrix `matrix`.
    """
    # Get the number of nodes and the edges for Gr
    B = m @ m.T  # The row graph of `m` is the undirected graph of `B = M @ M.T`
    nb_nodes = len(B)
    edges_list = np.argwhere(
        np.triu(B, k=1)
    )  # Get the indices of the non-zero values of the upper triangle of `B` without the diagonal

    # Generate the row-graph using `networkx`
    Gr = nx.Graph()
    nodes = range(nb_nodes)
    Gr.add_nodes_from(nodes, tag="inactive")  # Start with all the rows being 'inactive'
    Gr.add_edges_from(edges_list)
    return Gr


def get_nodes_s_and_e(Gr: nx.Graph) -> tuple[int, int]:
    """
    Purpose
    -------
    Get the starting node `s` and the ending node `e` out of the row graph `Gr` of
    the initial given matrix. They are normally obtained from the pseudo-diameter
    of `Gr` but right now, I use the diameter itself to find them.

    * If I find an algorithm to calculate the pseudo-diameter efficiently, it would
    be a good thing because at some point, calculating the diameter will probably
    become pretty costly...

    Parameters
    ----------
    `Gr` : networkx.Graph()
        The row graph of the initial given matrix.

    Returns
    -------
    `connected_pairs[random_id]` : tuple
        The starting (s) and the ending (e) nodes for the MSRO algorithm, chosen at
        random if there are more than one pairs that are separated by the same
        distance (the diameter of the graph)
    """
    # Calculate the diameter of the graph
    diameter = nx.algorithms.distance_measures.diameter(Gr)

    # Find the nodes that are at a distance equal to the diameter of the graph
    diameter_nodes = nx.algorithms.distance_measures.periphery(Gr)

    # Check for connections with distance d
    connected_pairs = []
    for n1 in diameter_nodes:
        for n2 in diameter_nodes:
            if (
                n1 < n2
                and nx.shortest_path_length(Gr, source=n1, target=n2) == diameter
            ):
                connected_pairs.append((n1, n2))

    # Get any of the 2 nodes that are separated by a distance = diameter
    # Maybe just return the first pair found in the nested `for` loops ?
    random_id = np.random.randint(len(connected_pairs))
    return connected_pairs[random_id]


def get_rcgain_i(m: np.ndarray, nb_already_assembled_rows: int, row_id: int) -> int:
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

    `nb_already_permuted_rows` : int
        The number of rows that have already been permuted in the matrix.

    `rows_to_permute` : list
        The 2 rows to permute in `m`.

    Returns
    -------
    `rcgain` : int
        The increase of the front row size and the column row size from the
        permutation of the 2 given rows of the matrix.

    `nold` : int
        The number of variables in row `i` that are candidates for elimination.
    """
    # Assemble the given row
    matrix = np.array(m, copy=True)
    row = matrix[row_id]
    matrix = np.delete(matrix, row_id, axis=0)
    matrix = np.insert(matrix, nb_already_assembled_rows, row, axis=0)

    # Find the values of `s` and `newc
    s = 0
    newc = 0
    for col_id in range(matrix.shape[1]):
        column = matrix[:, col_id]
        not_in_front = column[nb_already_assembled_rows + 1 :]
        if not np.any(not_in_front):
            s += 1
        in_front = column[: nb_already_assembled_rows + 1]
        if np.all(in_front[:-1] == 0) and in_front[-1] == 1:
            newc += 1

    # Evaluate rcgain for this row permutation (using the equation given in the paper)
    rgain = 1 - s
    cgain = newc - s
    rcgain = rgain + cgain
    nold = s
    return rcgain, nold


def initialize_P(m: np.ndarray, Gr: nx.Graph, e: int) -> list:
    """
    Purpose
    -------
    Evaluate the initial values of P (priority function) for all the rows of `m`.

    Parameters
    ----------
    `m` : np.ndarray
        The matrix on which we want to apply the MSRO algorithm.

    `Gr` : nx.Graph
        The row-graph of this matrix.

    `e` : int
        The ending node of the row-graph found from its diameter.

    Returns
    -------
    `P` : list
        The initial priority function, linked with `matrix` and `Gr`.

    `v` : float
        The normalizing factor in the expression of the priority function.
    """
    W1 = 2
    W2 = 1
    W3 = 0.2
    P = []
    distances = {}

    # Save the distances between each nodes and the target node in a dictionnary
    for i in range(m.shape[0]):
        distances[i] = nx.shortest_path_length(Gr, source=i, target=e)

    # Get the normalizing factor
    v = 1 / np.sqrt(np.sum(value**2 for value in distances.values()))

    # Evaluate the initial priority function
    for i in range(m.shape[0]):
        rcgain_i, nold_i = get_rcgain_i(m, 0, i)
        P_i = W1 * rcgain_i + W2 * v * distances[i] - W3 * nold_i
        P.append(P_i)
    return P, v


def get_row_to_assemble(
    Gr: nx.Graph,
    P: list,
    s: int,
) -> tuple[int, int]:
    """
    Purpose
    -------
    Get the row to assemble in the matrix by looking at the maximum value of the
    active rows in `P`.

    Parameters
    ----------
    `Gr` : nx.Graph
        The current row graph of the matrix for which we want to apply the
        MSRO algorithm.

    `P` : np.ndarray
        The current priority function of the algoritmh.

    `s` : int
        The starting row, so the first one to permute. In the graph, it is the
        first node in `(s,e)` that are at a distance `D` (the diameter of `Gr`).

    Returns
    -------
    `assembling_idx` : int
        The index where the chosen row will be placed in the matrix.

    `row_to_assemble_idx` : int
        The chosen row index to assemble in `m`.
    """
    nb_assembled_rows = 0
    for node in Gr.nodes():
        if Gr.nodes[node]["tag"] == "assembled":
            nb_assembled_rows += 1
    active_P_ids = [node for node in Gr.nodes() if Gr.nodes[node]["tag"] == "active"]
    if len(active_P_ids) > 0 and nb_assembled_rows > 0:
        P_copy = np.array(P, copy=True)
        not_active_nodes_ids = np.setdiff1d(
            np.arange(len(P_copy)), np.array(active_P_ids)
        )
        P_copy[not_active_nodes_ids] = np.NINF
        max_active_P_id = np.argmax(P_copy)
        assembling_idx = nb_assembled_rows
        row_to_assemble_idx = int(max_active_P_id)
    else:
        assembling_idx = 0
        row_to_assemble_idx = s
    return assembling_idx, row_to_assemble_idx


# Find a way to not have to update the row graph and work with less objects to upgrade
def update_Gr(Gr: nx.Graph, row_to_assemble_idx: int, assembling_idx: int) -> nx.Graph:
    """
    Purpose
    -------
    Relabel the nodes in the row graph so that it follows the modifications applied
    to the matrix `m` on which we apply the MSRO algorithm.

    Paramteres
    ----------
    `Gr` : nx.Graph
        The previous row graph of the matrix.

    `row_to_assemble_idx` : int
        The chosen row index to assemble next in the matrix.

    `assembling_idx` : int
        The row index where the chosen row will be placed.

    Returns
    -------
    `Gr` : nx.Graph
        The row graph with its nodes relabeled. All the nodes have their respetive
        tags also updated and their tags follow the node through the relabeling.
    """
    # Build the mapping for the relabeling of the nodes in Gr since the rows got moved in the matrix
    mapping = {row_to_assemble_idx: assembling_idx}
    if row_to_assemble_idx > assembling_idx:
        for index in range(assembling_idx, row_to_assemble_idx):
            mapping[index] = index + 1
    else:
        mapping[assembling_idx] = row_to_assemble_idx

    # Associate each neighbor nodes to the assembled row the state of 'active'
    for node in Gr.nodes():
        if Gr.nodes[node]["tag"] != "assembled":  # assembled rows don't change anymore
            if node in Gr.neighbors(row_to_assemble_idx):
                Gr.nodes[node]["tag"] = "active"
    # Relabel the nodes that were influenced by the assembling of the chosen row
    Gr = nx.relabel.relabel_nodes(Gr, mapping, copy=True)

    # Make the tags follow the corresponding nodes (it doesn't do that automatically)
    tags = []
    for node in Gr.nodes():
        if node in mapping.items():
            original_node = mapping[node]  # Get the original node based on the mapping
            tag = Gr.nodes[original_node]["tag"]  # Retrieve the original tag
            tags.append(tag)
    for node, tag in zip(Gr.nodes(), tags):
        Gr.nodes[node]["tag"] = tag

    # Specify the fact that the assembled row has been assembled
    Gr.nodes[assembling_idx]["tag"] = "assembled"
    return Gr


def update_m(
    m: np.ndarray, row_to_assemble_idx: int, assembling_idx: int
) -> np.ndarray:
    """
    Purpose
    -------
    Update the given matrix by assembling the chosen row at the chosen row index.

    Parameters
    ----------
    `m` : np.ndarray
        The given matrix on which we want to replace the chosen row at the given
        row index.

    `row_to_assemble_idx` : int
        The index of the chosen row to assemble.

    `assembling_idx` : int
        The row index to which the chosen row will go.

    Returns
    -------
    `m` : np.ndarray
        The updated matrix.
    """
    row = m[row_to_assemble_idx]
    m = np.delete(m, row_to_assemble_idx, axis=0)
    m = np.insert(m, assembling_idx, row, axis=0)
    return m


def update_P(
    m: np.ndarray,
    P: list,
    v: float,
    Gr: nx.Graph,
    row_to_assemble_idx: int,
    assembling_idx: int,
    e: int,
) -> list:
    """
    Purpose
    -------
    Update the values in the priority function. Only the values that are linked
    with active rows are updated, because the other ones are not taken into
    consideration if they are inactive.

    Parameters
    ----------
    `m` : np.ndarray
        The given matrix on which we want to apply te MSRO algorithm.

    `P` : list
        The previous priority function values.

    `v` : float
        The normalizing factor needed to evaluate the priority function.

    `Gr` : nx.Graph
        The row graph that has been properly updated according to the
        modifications done to the given matrix.

    `row_to_assemble_idx` : int
        The index of the chosen row to be assembled.

    `assembling_idx` : int
        The row index at which the chosen row will be placed.

    `e` : int
        The ending node in the row graph (the target node).

    Returns
    -------
    `P` : list
        The updated priority function values, where the assembled row is now
        linked with the value `np.NINF` and only the `active` rows values
        have been updated.
    """
    P = np.delete(P, row_to_assemble_idx)
    P = np.insert(P, assembling_idx, np.NINF)
    W1 = 2
    W2 = 1
    W3 = 0.2
    for i in Gr.nodes():
        if Gr.nodes[i]["tag"] == "active":
            rcgain_i, nold_i = get_rcgain_i(m, assembling_idx + 1, i)
            d_ie = nx.shortest_path_length(Gr, source=i, target=e)
            P_i = W1 * rcgain_i + W2 * v * d_ie - W3 * nold_i
            P[i] = P_i
    return P


def update_data(
    m: np.ndarray,
    P: list,
    v: float,
    Gr: nx.Graph,
    rows_to_switch: list,
    e: int,
) -> list:
    """
    Purpose
    -------
    This function lets us permute the desired rows in `m`, update the values in
    `P` and flag nodes in `Gr` that become inactive by specifying a value of
    `np.NINF` in `P` once the row has been permuted and it relabels the nodes
    in the row-graph `Gr` according to the row permutation.

    Parameters
    ----------
    `m` : np.ndarray
        The matrix on which to apply the MSRO algorithm.

    `P` : list
        The priority function applied in the MSRO algorithm.

    `v` : float
        The normalizing factor needed to evaluate the priority function.

    `Gr` : nx.Graph
        The row-graph of the matrix `m`.

    `rows_to_switch` : list
        The index in which the row to assemble will go and the row to assemble in `m`
        according to the values in `P`. Get them by using `get_row_to_assemble()`.

    `e` : int
        The 'ending' node found from the pseudo-diameter of the initial row-graph `Gr`.

    Returns
    -------
    `m` ; np.ndarray
        The matrix on which we want to apply the MSRO algorithm after 1 row permutation.

    `P` : np.array
        The updated priority function.

    `Gr` : nx.Graph
        The row-graph after relabeling the nodes because of the row permutation.

    `e` : int
        The updated ending node label. It needs to be changed since if the row id `e`
        has been selected by the priority function, then this row id will become
        `assembling_idx`, which is important to keep in memory.
    """
    assembling_idx, row_to_assemble_idx = rows_to_switch
    # Update the row graph (relabel the selected nodes and specify which nodes became active)
    Gr = update_Gr(Gr, row_to_assemble_idx, assembling_idx)
    if row_to_assemble_idx == e:
        e = assembling_idx

    # Update `m` (assemble the selected row in `m`)
    m = update_m(m, row_to_assemble_idx, assembling_idx)

    # Update the values in P
    P = update_P(m, P, v, Gr, row_to_assemble_idx, assembling_idx, e)

    return m, P, Gr, e


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
) -> np.ndarray:
    """
    Purpose
    -------
    Apply the [MSRO](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1506(199904/05)6:3%3C189::AID-NLA160%3E3.0.CO;2-C)
    algorithm on the given matrix `m` in order to minimize the row front size and
    the front column size at the same time. This algorithm can permute the rows or
    the columns of `m`. The weights `W1` and `W2` are set to (`W1`,`W2`) = (2,1) for
    the priority function `P`.

    Parameters
    ----------
    `m` : np.ndarray
        The matrix on which we want to apply the MSRO algorithm.

    `perm` : str
        Whether to permute the rows or the columns of the given matrix. Set to
        `columns` by default for the sites relabeling in an MPS-MPO problem.

    `show_rowGraph` : bool
        Whether to show (True) the row graph or not (False).

    Returns
    -------
    * Should return the new row/column order instead of the matrix

    `m` : np.ndarray
        The modified version of `m` after the application of the algorithm.

    Example
    -------
    >>> import numpy as np
    >>> from MatRexAlgs.msro import msro
    >>> m = np.random.randint(2, (6,6))
    >>> optimized_m = msro(m, perm = "columns")
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
    P, v = initialize_P(matrix, Gr, e)

    # Update those variables until P contains only -inf
    while not np.all(P == np.NINF):
        assembling_idx, row_to_assemble_idx = get_row_to_assemble(Gr, P, s)
        rows_to_switch = [assembling_idx, row_to_assemble_idx]
        matrix, P, Gr, e = update_data(matrix, P, v, Gr, rows_to_switch, e)

    # Return the matrix in the same format it was taken as an input
    if perm == "columns":
        return matrix.T
    elif perm == "rows":
        return matrix
