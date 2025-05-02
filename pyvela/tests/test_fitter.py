import os
from pyvela.fitter import VelaFitter
from pint.models import get_model_and_toas


def test_fitter():
    datadir = f"{os.path.dirname(os.path.realpath(__file__))}/datafiles"
    dataset = "NGC6440E"
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"
    m, t = get_model_and_toas(parfile, timfile, planets=True)
    ftr = VelaFitter(t, m)
    ftr.fit_toas()
    assert all([ftr.model[pname].uncertainty_value > 0 for pname in m.free_params])
