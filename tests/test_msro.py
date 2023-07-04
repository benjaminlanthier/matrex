"""
File for testing the functions in `msro.py` to assure that it works as intended.
"""
import sys, os

sys.path.append(os.path.realpath(os.path.dirname(__file__) + "/.."))

import numpy as np
import pytest

from matrex.msro import (  # type: ignore
    find_active_rows,
    change_in_frontsize,
    find_front_columns,
    calculate_ordering,
    msro,
)


def test_find_active_rows():
    """Test for the `matrex.msro.test_find_active_rows()` function."""
    reordered = np.array([3, 1])
    matrix = np.array(
        [
            [1, 0, 1, 1, 0, 0],
            [0, 1, 0, 1, 1, 0],
            [1, 0, 1, 1, 0, 1],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 1, 1],
            [0, 0, 0, 0, 0, 1],
        ]
    )
    row_graph = matrix @ matrix.T
    active_rows = find_active_rows(reordered, row_graph)
    assert np.all(active_rows == np.array([0, 2, 4]))


def test_change_in_frontsize():
    """Test for the `matrex.msro.change_in_frontsize()` function."""
    # Test for this matrix with 2 rows in the front
    # matrix = np.array([[0, 1, 1, 0, 0],
    #                    [1, 0, 0, 1, 0],
    #                    [1, 0, 0, 1, 1],
    #                    [1, 1, 1, 1, 0],
    #                    [0, 1, 0, 1, 0],
    #                    [1, 0, 1, 1, 1]])
    front = np.array([[0, 1, 1, 0, 0], [1, 0, 0, 1, 0]])
    candidate_rows = np.array(
        [[1, 0, 0, 1, 1], [1, 1, 1, 1, 0], [0, 1, 0, 1, 0], [1, 0, 1, 1, 1]]
    )
    remaining_rows = np.array(
        [
            [[1, 1, 1, 1, 0], [0, 1, 0, 1, 0], [1, 0, 1, 1, 1]],
            [[1, 0, 0, 1, 1], [0, 1, 0, 1, 0], [1, 0, 1, 1, 1]],
            [[1, 0, 0, 1, 1], [1, 1, 1, 1, 0], [1, 0, 1, 1, 1]],
            [[1, 0, 0, 1, 1], [1, 1, 1, 1, 0], [0, 1, 0, 1, 0]],
        ]
    )
    rcgain = change_in_frontsize(front, candidate_rows, remaining_rows)
    assert np.all(rcgain == np.array([2, 1, 1, 2]))


def test_find_front_columns():
    """Test for the `matrex.msro.find_front_columns()` function."""
    # TODO
    # Test for this matrix with 2 rows in the front
    # matrix = np.array([[0, 1, 1, 0, 0],
    #                    [1, 0, 0, 1, 0],
    #                    [1, 0, 0, 1, 1],
    #                    [1, 1, 1, 1, 0],
    #                    [0, 1, 0, 1, 0],
    #                    [1, 0, 1, 1, 1]])
    front = np.array([[0, 1, 1, 0, 0], [1, 0, 0, 1, 0]])
    candidate_rows = np.array(
        [[1, 0, 0, 1, 1], [1, 1, 1, 1, 0], [0, 1, 0, 1, 0], [1, 0, 1, 1, 1]]
    )
    front_columns = find_front_columns(front, candidate_rows)
    assert np.all(front_columns == np.array([2, 4, 2, 3]))


@pytest.mark.parametrize("m, n", [(10, 10), (20, 50), (50, 20), (100, 100)])
def test_msro_result1(m, n):
    """Test to validate the `matrex.msro.msro()` function's output."""
    matrix = np.zeros((m, n), dtype=int)
    for i in range(m):
        indices = np.random.choice(a=range(n), size=3, replace=False)
        matrix[i, indices] = 1
    new_rows_order = msro(input_matrix=matrix)
    assert len(new_rows_order) == len(set(new_rows_order))


def test_msro_result2():
    """Test to validate the `matrex.msro.msro()` function's output."""
    matrix = np.array(
        [
            [1, 0, 1, 1, 0, 0],
            [0, 1, 0, 1, 1, 0],
            [1, 0, 1, 1, 0, 1],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 1, 1],
            [0, 0, 0, 0, 0, 1],
        ]
    )
    new_rows_order = msro(input_matrix=matrix)
    goal_rows_order = [3, 1, 4, 5, 2, 0]
    assert new_rows_order == goal_rows_order
