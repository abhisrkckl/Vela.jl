from pint2vela.ecorr import ecorr_sort
from pint.models import get_model_and_toas
from pint.config import examplefile

import numpy as np
import pytest


@pytest.fixture
def data_B1855p09():
    parfile = examplefile("B1855+09_NANOGrav_9yv1.gls.par")
    timfile = examplefile("B1855+09_NANOGrav_9yv1.tim")
    return get_model_and_toas(parfile, timfile)


def test_ecorr_sort(data_B1855p09):
    m, t = data_B1855p09

    t_sorted, ranges, indices = ecorr_sort(m, t)

    assert len(t_sorted) == len(t)
    assert np.all(np.diff(t_sorted.table["index"]) == 1)

    assert len(ranges) == len(indices)

    for i in range(1, len(ranges)):
        assert ranges[i][0] == ranges[i - 1][1] + 1

    assert (
        len(np.unique(m.components["EcorrNoise"].get_noise_weights(t)))
        == len(np.unique(indices)) - 1
    )
