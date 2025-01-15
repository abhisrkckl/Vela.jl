from io import StringIO
import os

import emcee
import numpy as np
import pytest

from pyvela.spnta import SPNTA
from pyvela import (
    pyvela_compare_script,
    pyvela_jlso_script,
    pyvela_script,
    pyvela_plot_script,
)
from pint.models import get_model_and_toas

prior_str = """
    {
        "EFAC": {
            "distribution": "Uniform",
            "args": [0.5, 1.5]
        },
        "EQUAD": {
            "distribution": "Normal",
            "args": [0.0, 1.0],
            "lower": 0.0
        }
    }
"""


def test_analysis_NGC6440E_emcee():
    datadir = os.path.dirname(os.path.realpath(__file__)) + "/datafiles"
    dataset = "NGC6440E"
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"

    mp, tp = get_model_and_toas(parfile, timfile, planets=True)

    spnta = SPNTA.from_pint(
        mp,
        tp,
        custom_priors=StringIO(prior_str),
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


@pytest.mark.parametrize("dataset", ["NGC6440E", "sim_sw.wb"])
def test_script(dataset):
    datadir = os.path.dirname(os.path.realpath(__file__)) + "/datafiles"
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"
    outdir = f"_{dataset}_out"

    prior_file = "__prior.json"
    with open(prior_file, "w") as pf:
        print(prior_str, file=pf)

    args = f"{parfile} {timfile} -P {prior_file} -T {parfile} -o {outdir}".split()

    pyvela_script.main(args)

    assert os.path.isdir(outdir)
    assert os.path.isfile(f"{outdir}/summary.json")
    assert os.path.isfile(f"{outdir}/{prior_file}")
    assert os.path.isfile(f"{outdir}/param_names.txt")
    assert os.path.isfile(f"{outdir}/param_units.txt")
    assert os.path.isfile(f"{outdir}/param_scale_factors.txt")
    assert os.path.isfile(f"{outdir}/samples_raw.npy")
    assert os.path.isfile(f"{outdir}/samples.npy")

    pyvela_plot_script.main([f"{outdir}/"])


@pytest.mark.parametrize("dataset", ["NGC6440E", "sim_sw.wb"])
def test_compare_script(dataset):
    datadir = os.path.dirname(os.path.realpath(__file__)) + "/datafiles"
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"

    prior_file = "__prior.json"
    with open(prior_file, "w") as pf:
        print(prior_str, file=pf)

    args = f"{parfile} {timfile} -P {prior_file}".split()

    pyvela_compare_script.main(args)


@pytest.mark.parametrize("dataset", ["NGC6440E"])
def test_jlso_script(dataset):
    datadir = os.path.dirname(os.path.realpath(__file__)) + "/datafiles"
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"

    prior_file = "__prior.json"
    with open(prior_file, "w") as pf:
        print(prior_str, file=pf)

    outfile = f"__{dataset}.jlso"

    args = f"{parfile} {timfile} -P {prior_file} -o {outfile}".split()

    pyvela_jlso_script.main(args)

    assert os.path.isfile(outfile)
