import numpy as np
import pytest

from pint2vela import read_model_and_toas
from pint2vela.vela import vl


@pytest.fixture
def data_sim_jump():
    mv, tv = read_model_and_toas("datafiles/sim_sw.par", "datafiles/sim_sw.tim")
    params = vl.read_param_values_to_vector(
        mv.param_handler, mv.param_handler._default_params_tuple
    )
    return mv, tv, params


def test_data(data_sim_jump):
    mv, tv, params = data_sim_jump
    assert len(tv) == 1000
    assert len(mv.components) == 5
    assert set(vl.get_free_param_names(mv.param_handler)) == {
        "RAJ",
        "DECJ",
        "DM",
        "PHOFF",
        "F0",
        "F1",
        "NE_SW"
    }
    assert len(params) == 8


def test_chi2(data_sim_jump):
    mv, tv, params = data_sim_jump
    calc_chi2 = vl.get_chi2_func(mv, tv)
    # assert calc_chi2(params) / len(tv) < 1.1


def test_likelihood(data_sim_jump):
    mv, tv, params = data_sim_jump
    calc_lnlike = vl.get_lnlike_func(mv, tv)
    assert np.isfinite(calc_lnlike(params))


def test_prior(data_sim_jump):
    mv, _, params = data_sim_jump
    calc_lnprior = vl.get_lnprior_func(mv)
    assert np.isfinite(calc_lnprior(mv.param_handler._default_params_tuple))
    assert calc_lnprior(params) == calc_lnprior(mv.param_handler._default_params_tuple)
