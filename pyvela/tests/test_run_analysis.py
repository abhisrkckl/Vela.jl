import os

import emcee
import numpy as np

from pyvela.spnta import SPNTA


def test_analysis_NGC6440E_emcee():
    datadir = os.path.dirname(os.path.realpath(__file__)) + "/datafiles"
    dataset = "NGC6440E"
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"

    spnta = SPNTA(
        parfile,
        timfile,
    )

    nwalkers = 3 * spnta.ndim

    p0 = np.array(
        [spnta.prior_transform(cube) for cube in np.random.rand(nwalkers, spnta.ndim)]
    )

    sampler = emcee.EnsembleSampler(nwalkers, spnta.ndim, spnta.lnpost)
    sampler.run_mcmc(p0, 1500)

    samples_raw = sampler.get_chain(flat=True, discard=500, thin=50)
    samples = spnta.rescale_samples(samples_raw)

    pint_model_final = spnta.update_pint_model(samples)

    assert set(pint_model_final.free_params) == set(spnta.model_pint.free_params)
    assert all(
        pint_model_final[pname].uncertainty is not None
        and pint_model_final[pname].uncertainty_value != 0
        for pname in pint_model_final.free_params
    )
