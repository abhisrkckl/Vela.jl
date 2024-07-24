import numpy as np
from pint2vela import read_model_and_toas, vl
import pytest

@pytest.fixture
def data_pure_rotator():
    mv, tv = read_model_and_toas("datafiles/sim_jump.par", "datafiles/sim_jump.tim")
    params = vl.read_param_values_to_vector(mv.param_handler, mv.param_handler._default_params_tuple)
    return mv, tv, params

def test_data(data_pure_rotator):
    mv, tv, params = data_pure_rotator
    assert len(tv) == 1000
    assert len(mv.components) == 5
    assert set(vl.get_free_param_names(mv.param_handler)) == {"RAJ", "DECJ", "DM", "PHOFF", "F0", "F1", "JUMP1", "JUMP2"}
    assert len(params) == 8

def test_chi2(data_pure_rotator):
    mv, tv, params = data_pure_rotator
    calc_chi2 = vl.get_chi2_func(mv, tv)
    assert calc_chi2(params) / len(tv) < 1.1

def test_likelihood(data_pure_rotator):
    mv, tv, params = data_pure_rotator
    calc_lnlike = vl.get_lnlike_func(mv, tv)
    assert np.isfinite(calc_lnlike(params))