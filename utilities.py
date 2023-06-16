"""
Purpose: A set of efficient, matrix-based tools for calculating the reordering of a
         matrix's rows. As much as possible, I use matrix-based computations instead
         of loops.
Date created: 2023-06-14
"""

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
  frontColumns = np.any(a = front, axis = 0) # Boolean mask for which indices are present
  newColumns = np.sum(np.logical_and(~frontColumns, row))
  fullySummed = np.sum(np.logical_and(~np.any(a = remaining, axis = 0), row))
  return 1 + newColumns - 2 * fullySummed
  
def pseudodiameter(graph):
  """
    - Purpose: Find nodes i and j which have a distance close to the diameter.
               Based on Mathematica's approach:
               https://reference.wolfram.com/language/GraphUtilities/ref/PseudoDiameter.html
    - Input:
      - graph (NetworkX object): The row graph of the matrix.
    - Output:
      - (i, j) (tuple): The node indices that are separated by the pseudodiameter.
  """
  nodes = graph.nodes()
  start = np.random.choice(nodes)
  finalLoop = False
  distance = 0
  while True:
    distances = nx.algorithms.single_source_shortest_path_length(G = graph, source = start)
    farthestNode = max(distances, key = distances.get)
    newDistance = distances[farthestNode]
    if finalLoop:
      break
    if newDistance == distance:
      finalDegrees = dict(graph.degree(nbunch = distances.keys()))
      start = min(finalDegrees, key = finalDegrees.get)
      finalLoop = True
    else:
      start = farthestNode
      distance = newDistance

  return (start, farthestNode)

def findActiveRows(reordered, rowGraph):
  """
    - Purpose: Find row indices which correspond to active rows.
    - Input:
      - reorderd (list of integers): The current nodes that are part of the front.
      - rowGraph (array): The row graph of the matrix.
    - Output:
      - active (list of integers): The active node indices.
  """
  # Use the row graph to find indices of active nodes that aren't already part of the front.

  candidates = np.nonzero(np.any(rowGraph[reordered], axis = 0))[0]
  active = np.setdiff1d(ar1 = candidates, ar2 = reordered)
  return active

def findFrontColumns(matrix, row):
  """
    - Purpose: For a given row, calculate the number of columns that
               have nonzero elements and are already in the front.
    - Input:
      - matrix (array): The current matrix that describes the front.
      - row (array): The particular row to consider.
    - Output:
      - (integer): The number of columns already in the front.
  """
  currentFrontColumns = np.any(a = matrix, axis = 0) # Boolean mask
  consider = np.any(row)
  return np.sum(np.logical_and(~currentFrontColumns, consider))

def calculateOrdering(graph, matrix, W1 = 2, W2 = 1, W3 = 0.2, verbose = False):
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
    for node in active:
      order.append(node)
      unordered = np.setdiff1d(ar1 = allNodes, ar2 = order)
      order.pop()
      row, remaining = matrix[node], matrix[unordered]
      frontColumns = findFrontColumns(matrix = front, row = row)
      deltaFront = change_in_frontsize(front = front, row = row, remaining = remaining)
      score = -W1 * deltaFront + W2 * distances[node] - W3 * frontColumns
      priorities[node] = score
      if verbose:
        print("Node: ", node)
        print("Distance to target: ", distances[node])
        print("Columns already in the front: ", findFrontColumns(matrix = front, row = row))
        print("Order    : ", order)
        print("Unordered: ", unordered)
    # Choose highest priority
    selection = max(priorities, key = priorities.get)
    order.append(selection)
    if verbose:
      print("Start, target: ", start, target)
      print("Active: ", active)
      print("Priorities: ")
      print(priorities)
      print("Selection: ", selection)
      print()
  assert len(np.unique(order)) == len(allNodes)
  return order

def MSRO(input_matrix, W1 = 2, W2 = 1, W3 = 0.2, verbose = False):
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
    order = calculateOrdering(graph = graph, matrix = matrix, W1 = W1, W2 = W2, W3 = W3, verbose = verbose)
    totalReordering.extend(order)
  return totalReordering