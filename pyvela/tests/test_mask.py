import numpy as np
import pytest

from pyvela.model import get_exclusive_mask, is_exclusive_mask


def test_non_exclusive_mask():
    mask1 = np.array([[1, 1, 0], [0, 1, 0]])
    assert not is_exclusive_mask(mask1)
    with pytest.raises(ValueError):
        get_exclusive_mask(mask1)


def test_exclusive_mask():
    mask1 = np.array([[1, 0, 0], [0, 1, 0]])
    assert is_exclusive_mask(mask1)
    exmask = get_exclusive_mask(mask1)
    assert all(exmask == np.array([1, 2, 0]))
