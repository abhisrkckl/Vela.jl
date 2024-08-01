import numpy as np
import pytest

from pint2vela import read_model_and_toas
from pint2vela.vela import vl


@pytest.fixture
def data_sim1():
    mv, tv = read_model_and_toas("datafiles/sim1.par", "datafiles/sim1.tim")
    params = vl.read_param_values_to_vector(
        mv.param_handler, mv.param_handler._default_params_tuple
    )
    return mv, tv, params


def test_data(data_sim1):
    mv, tv, params = data_sim1
    assert len(tv) == 2000
    assert len(mv.components) == 5
    assert set(vl.get_free_param_names(mv.param_handler)) == {
        "RAJ",
        "DECJ",
        "DM",
        "PHOFF",
        "F0",
        "F1",
        "EFAC1",
        "EQUAD1",
    }
    assert len(params) == 8


def test_chi2(data_sim1):
    mv, tv, params = data_sim1
    calc_chi2 = vl.get_chi2_func(mv, tv)
    assert calc_chi2(params) / len(tv) < 1.1


def test_likelihood(data_sim1):
    mv, tv, params = data_sim1
    calc_lnlike = vl.get_lnlike_func(mv, tv)
    assert np.isfinite(calc_lnlike(params))


def test_prior(data_sim1):
    mv, _, params = data_sim1
    calc_lnprior = vl.get_lnprior_func(mv)
    assert np.isfinite(calc_lnprior(mv.param_handler._default_params_tuple))
    assert calc_lnprior(params) == calc_lnprior(mv.param_handler._default_params_tuple)
