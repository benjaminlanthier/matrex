"""
Purpose: Demonstrate an example of the reordering algorithm for a random matrix.

Date created: 2023-06-16
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import networkx as nx
from utilities import MSRO
from vectorizedUtilities import vectorizedMSRO

def getFirstOrdering(matrix):
    """
        - Purpose: Find the first 
    """
    first = []
    for row in matrix:
        first.append(np.nonzero(row)[0][0])
    return first

def getLastOrdering(matrix):
    last = []
    for row in matrix:
        last.append(np.nonzero(row)[0][-1])
    return last

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
    N = 300
    M = 300
    matrix = np.zeros((M,N), dtype = int)
    for i in range(M):
        indices = np.random.choice(a = range(N), size = 3, replace = False)
        matrix[i, indices] = 1
    
    W1, W2, W3 = 2, 1, 0.2
    algorithmTime = time.time()
    #order = MSRO(input_matrix = matrix.T, W1 = W1, W2 = W2, W3 = W3)
    order = vectorizedMSRO(input_matrix = matrix.T, W1 = W1, W2 = W2, W3 = W3)
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

    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (10, 8))
    ax1.imshow(matrix, cmap = "binary")
    ax2.imshow(reordered, cmap = "binary")
    ax1.set_title(r"Original matrix ($W = {:.2f}$)".format(np.mean(widths)))
    ax2.set_title(r"Reduced matrix ($W = {:.2f}$)".format(np.mean(reorderedWidths)))
    plt.savefig("reorderedMatrixN{}.pdf".format(N), dpi = 300)
    #plt.show()