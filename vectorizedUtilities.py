"""
Purpose: A set of efficient, matrix-based tools for calculating the reordering of a
         matrix's rows. As much as possible, I use matrix-based computations instead
         of loops.
Date created: 2023-06-16
"""

import numpy as np
import networkx as nx
from utilities import pseudodiameter, findActiveRows
from operator import itemgetter

def change_in_frontsize_vectorized(front, rows, remainings):
  """
    - Purpose: Calculate the change in row and column frontsizes.
    - Input:
      - front (array): The front matrix before adding the current row.
      - rows (array of arrays): The new rows to consider.
      - remainings (array of arrays): The remaining rows.
    - Output:
      - (integer): The change in front sizes.
  """
  frontColumns = (~np.any(a = front, axis = 0)).astype(int)
  newColumns = rows.astype(int) @ frontColumns.T
  zerosInRemainings = (~np.any(a = remainings, axis = 1)).astype(int)
  fullySummed = np.sum(a = zerosInRemainings @ rows.astype(int).T, axis = 1)
  
  return 1 + newColumns - 2 * fullySummed

def findFrontColumnsVectorized(matrix, rows):
  """
    - Purpose: For a given row, calculate the number of columns that
               have nonzero elements and are already in the front.
    - Input:
      - matrix (array): The current matrix that describes the front.
      - rows (array of shape (number of matrices, len(matrix[0]))):
               A multidimensional array holding all of the considered rows.
    - Output:
      - frontColumns (array of integers): The number of columns already in the front.
  """
  currentFrontColumns = (~np.any(a = matrix, axis = 0)).astype(int)
  return rows.astype(int) @ currentFrontColumns.T

def calculateOrderingVectorized(graph, matrix, W1 = 2, W2 = 1, W3 = 0.2, verbose = False):
  """
    - Purpose: For a connected graph, calculate an ordering of the nodes.
    - Inputs:
      - graph (NetworkX object): The graph to reorder.
      - matrix (array): The matrix for the entire original graph.
      - W1, W2, W3 (positive integers): The weights for balancing the algorithm's
                                        considerations.
      - verbose (Boolean): If True, prints out information about the process.
    - Output:
      - order (list of integers): A permutation of the graph nodes.
  """
  rowGraph = matrix @ matrix.T
  start, target = pseudodiameter(graph = graph)
  distances = nx.algorithms.single_source_shortest_path_length(G = graph, source = target)
  if verbose:
    print("Distances: ", distances)
  allNodes = graph.nodes()
  order = [start]
  while len(order) < len(graph):
    # Find active rows
    active = findActiveRows(reordered = order, rowGraph = rowGraph)
    priorities = {}
    front = matrix[order]
    
    # Build the matrix of indices for each list of
    # ordered or unordered indices
    multiple = np.tile(A = order, reps = (len(active), 1))
    multiple = np.column_stack(tup = (multiple, active))
    
    unorderedIndices = np.setdiff1d(ar1 = allNodes, ar2 = order)
    
    # Still need to understand these sorting lines!
    removeMask = (unorderedIndices==active[:,None]).any(0).nonzero()[0]
    unorderedIndices = np.tile(A = unorderedIndices, reps = (len(active), 1))
    m, n = unorderedIndices.shape
    unorderedIndices = unorderedIndices[np.arange(n) != removeMask[:,None]].reshape(m,-1)
    
    # Build the rows and remaining matrices for each row of the orderedMask
    rows = matrix[active]
    remainings = matrix[unorderedIndices]
    
    # Calculate function values in matricized form
    values = findFrontColumnsVectorized(matrix = front, rows = rows)
    deltas = change_in_frontsize_vectorized(front = front, rows = rows, remainings = remainings)
    score = -W1 * deltas + W2 * itemgetter(*tuple(active))(distances) - W3 * values
    selectionIndex = np.argmax(a = score)
    selection = active[selectionIndex]
    order.append(selection)
    if verbose:
      print("Start, target: ", start, target)
      print("Active: ", active)
      print("Priorities: ")
      print(priorities)
      print("Selection: ", selection)
      print()
  assert len(np.unique(order)) == len(allNodes)
  #print("Ordering time: ", time.time() - orderTime)
  return order

def vectorizedMSRO(input_matrix, W1 = 2, W2 = 1, W3 = 0.2, verbose = False):
  """
    - Purpose: Reorder the rows of the input matrix to minimize its front.
    - Input:
      - input_matrix (array): The matrix to reorder.
      - W1, W2, W3 (positive integers): The weights for balancing the algorithm's
                                        considerations.
      - verbose (Boolean): If True, prints out information about the process.
    - Output:
      - totalReordering (1D array): The permutation of the rows.
  """

  # Only consider the nonzero pattern of elements
  matrix = input_matrix.astype(bool)
  
  # Compute initial distances between rows of the matrix.
  rowGraph = nx.from_numpy_array(matrix @ matrix.T)

  # Check if the graph is disconnected and break into components
  totalReordering = []
  subgraphs = [rowGraph.subgraph(component).copy() for component in nx.connected_components(rowGraph)]
  if verbose:
    print("Component sizes: ", [len(i) for i in subgraphs])
  for e, graph in enumerate(subgraphs):
    if verbose:
      print("Subgraph: ", e)
    order = calculateOrderingVectorized(graph = graph, matrix = matrix, W1 = W1, W2 = W2, W3 = W3, verbose = verbose)
    totalReordering.extend(order)
  return totalReordering