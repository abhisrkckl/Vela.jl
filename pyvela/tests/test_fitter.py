import os
import numpy as np
from pint.models import get_model_and_toas
from pyvela.fitter import VelaFitter


def test_vela_fitter():
    datadir = f"{os.path.dirname(os.path.realpath(__file__))}/datafiles"
    dataset = "NGC6440E"
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"

    mp, tp = get_model_and_toas(parfile, timfile)

    ftr = VelaFitter(tp, mp)

    ftr.fit_toas(mcmc=True)
    pvals_mcmc = np.array([ftr.model[par].value for par in ftr.model.free_params])
    perrs_mcmc = np.array(
        [ftr.model[par].uncertainty_value for par in ftr.model.free_params]
    )
    assert np.all(map(np.isfinite, pvals_mcmc))
    assert np.all(map(np.isfinite, perrs_mcmc))
    assert ftr.samples is not None

    ftr.fit_toas(mcmc=False)
    pvals_maxpost = np.array([ftr.model[par].value for par in ftr.model.free_params])
    perrs_maxpost = np.array(
        [ftr.model[par].uncertainty_value for par in ftr.model.free_params]
    )
    assert np.all(map(np.isfinite, pvals_maxpost))
    assert np.all(map(np.isfinite, perrs_maxpost))
    assert ftr.samples is None

    assert np.all((pvals_mcmc - pvals_maxpost) / perrs_mcmc < 1.5)
