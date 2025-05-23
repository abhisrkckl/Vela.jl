import os

import pytest
from pint.models import get_model_and_toas

from pyvela.fitter import VelaFitter


@pytest.mark.parametrize("dataset", ["NGC6440E", "sim_sw.wb"])
def test_fitter(dataset):
    datadir = f"{os.path.dirname(os.path.realpath(__file__))}/datafiles"
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"
    m, t = get_model_and_toas(parfile, timfile, planets=True)
    ftr = VelaFitter(t, m)
    ftr.fit_toas()
    assert all([ftr.model[pname].uncertainty_value > 0 for pname in m.free_params])
