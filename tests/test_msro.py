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
    # TODO
    pass


def test_change_in_frontsize():
    """Test for the `matrex.msro.change_in_frontsize()` function."""
    # TODO
    pass


def test_find_front_columns():
    """Test for the `matrex.msro.find_front_columns()` function."""
    # TODO

    pass


def test_calculate_ordering():
    """Test for the `matrex.msro.calculate_ordering()` function."""
    # TODO
    pass


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
