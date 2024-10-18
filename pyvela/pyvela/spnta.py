import numpy as np
from pint.binaryconvert import convert_binary
from pint.logging import setup as setup_log
from pint.models import get_model_and_toas

from .model import pint_model_to_vela
from .toas import pint_toas_to_vela
from .ecorr import ecorr_sort
from .vela import vl


def read_model_and_toas(
    parfile: str, timfile: str, cheat_prior_scale=20, custom_priors={}
):
    """Read a pair of par & tim files and create a `Vela.TimingModel` object and a
    Julia `Vector` of `TOA`s."""

    setup_log(level="WARNING")
    mp, tp = get_model_and_toas(
        parfile,
        timfile,
        planets=True,
        allow_tcb=True,
        allow_T2=True,
        add_tzr_to_model=True,
    )

    if "BinaryBT" in mp.components:
        mp = convert_binary(mp, "DD")

    if "EcorrNoise" in mp.components:
        assert not tp.is_wideband(), "ECORR is not supported for wideband data."
        tp, ecorr_toa_ranges, ecorr_indices = ecorr_sort(mp, tp)
    else:
        ecorr_toa_ranges, ecorr_indices = None, None

    model = pint_model_to_vela(
        mp,
        tp,
        cheat_prior_scale,
        custom_priors,
        ecorr_toa_ranges=ecorr_toa_ranges,
        ecorr_indices=ecorr_indices,
    )
    toas = pint_toas_to_vela(tp, float(mp.PEPOCH.value))

    return model, toas


class SPNTA:
    def __init__(
        self,
        parfile: str,
        timfile: str,
        cheat_prior_scale: float = 20,
        custom_priors: dict = {},
    ):
        self.parfile = parfile
        self.timfile = timfile

        model, toas = read_model_and_toas(
            parfile,
            timfile,
            cheat_prior_scale=cheat_prior_scale,
            custom_priors=custom_priors,
        )

        self._setup(model, toas)

    def _setup(self, model, toas):
        self.model = model
        self.toas = toas

        self.lnlike = vl.get_lnlike_func(self.model, self.toas)
        self.lnprior = vl.get_lnprior_func(self.model)
        self.prior_transform = vl.get_prior_transform_func(self.model)
        self.lnpost = vl.get_lnpost_func(self.model, self.toas, False)
        self.lnpost_vectorized = vl.get_lnpost_func(self.model, self.toas, True)

        self.param_names = list(vl.get_free_param_names(self.model))
        self.param_labels = list(vl.get_free_param_labels(self.model))
        self.scale_factors = np.array(vl.get_scale_factors(self.model))
        self.maxlike_params = np.array(vl.read_param_values_to_vector(self.model))
        self.ndim = len(self.param_names)

        self._check()

    def _check(self):
        cube = np.random.rand(self.ndim)
        sample = self.prior_transform(cube)
        lnpr = self.lnprior(sample)
        lnl = self.lnlike(sample)
        lnp = self.lnpost(sample)
        lnpv = self.lnpost_vectorized(np.array([sample]))
        assert all(np.isfinite([lnpr, lnl, lnp]))
        assert np.isclose(lnp, lnpv)

    def rescale_samples(self, samples: np.ndarray):
        return samples / self.scale_factors

    def save_jlso(self, filename: str):
        vl.save_pulsar_data(filename, self.model, self.toas)

    @classmethod
    def load_jlso(cls, filename: str):
        spnta = cls.__new__(cls)
        model, toas = vl.load_pulsar_data(filename)
        spnta.jlsofile = filename
        spnta._setup(model, toas)
        return spnta
