from copy import deepcopy
import json
from typing import IO

import numpy as np
import astropy.units as u
from pint.binaryconvert import convert_binary
from pint.logging import setup as setup_log
from pint.models import TimingModel, get_model_and_toas
from pint.toa import TOAs
from pint import dmu, DMconst

from .ecorr import ecorr_sort
from .model import pint_model_to_vela
from .priors import process_custom_priors
from .toas import day_to_s, pint_toas_to_vela
from .vela import jl, vl


def convert_model_and_toas(
    model: TimingModel, toas: TOAs, cheat_prior_scale=20, custom_priors={}
):
    """Read a pair of par & tim files and create a `Vela.TimingModel` object and a
    Julia `Vector` of `TOA`s."""

    if "BinaryBT" in model.components:
        model = convert_binary(model, "DD")

    if "EcorrNoise" in model.components:
        assert not toas.is_wideband(), "ECORR is not supported for wideband data."
        toas, ecorr_toa_ranges, ecorr_indices = ecorr_sort(model, toas)
    else:
        ecorr_toa_ranges, ecorr_indices = None, None

    model_v = pint_model_to_vela(
        model,
        toas,
        cheat_prior_scale,
        custom_priors,
        ecorr_toa_ranges=ecorr_toa_ranges,
        ecorr_indices=ecorr_indices,
    )
    toas_v = pint_toas_to_vela(toas, float(model.PEPOCH.value))

    return model_v, toas_v


class SPNTA:
    """
    A Python class that wraps the Julia objects representing the timing & noise model (`Vela.TimingModel`)
    and a collection of TOAs (`Vector{Vela.TOA}`). It provides various callable objects that can be used
    to run Bayesian inference.

    Parameters & Attributes
    -----------------------
    parfile: str
        Name of the par file from which the model was constructed.
    timfile: str
        Name of the tim file from which the TOAs were constructed.
    model_pint: pint.models.TimingModel
        A PINT `TimingModel` object
    model: Vela.TimingModel
        A Vela `TimingModel` object
    toas: Vector[Vela.TOA]
        A collection of Vela TOAs.
    lnlike: Callable
        Computes the log-likelihood
    lnprior: Callable
        Computes the log-prior
    prior_transform: Callable
        Computes the prior transform
    lnpost: Callable
        Computes the log-posterior
    lnpost_vectorized: Callable
        Computes the log-posterior for a collection of parameters
    param_names: List[str]
        Free parameter names
    param_labels: List[str]
        Free parameter labels (for plotting)
    scale_factors: numpy.ndarray
        Scale factors for converting between Vela's internal units and the commonly used units.
    maxlike_params: numpy.ndarray
        Free parameter values taken from the par file in internal units
    ndim: int
        Number of free parameters
    """

    def __init__(
        self,
        parfile: str,
        timfile: str,
        cheat_prior_scale: float = 20,
        custom_priors: str | IO | dict = {},
        pint_kwargs={},
    ):
        self.parfile = parfile
        self.timfile = timfile

        setup_log(level="WARNING")
        model_pint, toas_pint = get_model_and_toas(
            parfile,
            timfile,
            planets=True,
            allow_T2=True,
            allow_tcb=True,
            add_tzr_to_model=True,
            **pint_kwargs,
        )
        self.model_pint = model_pint

        # custom_priors_dict is in the "raw" format. The numbers may be
        # in "normal" units and have to be converted into internal units.
        if isinstance(custom_priors, dict):
            custom_priors_dict = custom_priors
        elif isinstance(custom_priors, str):
            with open(custom_priors) as custom_priors_file:
                custom_priors_dict = json.load(custom_priors_file)
        else:
            custom_priors_dict = json.load(custom_priors)

        custom_priors = process_custom_priors(custom_priors_dict, model_pint)

        setup_log(level="WARNING")
        model, toas = convert_model_and_toas(
            model_pint,
            toas_pint,
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

    def rescale_samples(self, samples: np.ndarray) -> np.ndarray:
        """Rescale the samples from Vela's internal units to common units"""
        return samples / self.scale_factors

    def save_jlso(self, filename: str) -> None:
        """Write the model and TOAs as a JLSO file"""
        vl.save_pulsar_data(filename, self.model, self.toas)

    def is_wideband(self) -> bool:
        """Whether the TOAs are wideband."""
        return jl.isa(self.toas[0], vl.WidebandTOA)

    def get_mjds(self) -> np.ndarray:
        """Get the MJDs of each TOA."""
        if self.is_wideband():
            return np.array(
                [jl.Float64(vl.value(wtoa.toa.value)) / day_to_s for wtoa in self.toas]
            )
        else:
            return np.array(
                [jl.Float64(vl.value(toa.value)) / day_to_s for toa in self.toas]
            )

    def time_residuals(self, params: np.ndarray) -> np.ndarray:
        """Get the timing residuals (s) for a given set of parameters."""
        params = vl.read_params(self.model, params)
        return np.array(
            list(map(vl.value, vl.form_residuals(self.model, self.toas, params)))
            if not self.is_wideband()
            else [
                vl.value(wr[0])
                for wr in vl.form_residuals(self.model, self.toas, params)
            ]
        )

    def scaled_toa_unceritainties(self, params: np.ndarray) -> np.ndarray:
        """Get the scaled TOA uncertainties (s) for a given set of parameters."""
        params = vl.read_params(self.model, params)

        ctoas = [vl.correct_toa(self.model, tvi, params) for tvi in self.toas]

        return np.sqrt(
            [
                vl.value(
                    vl.scaled_toa_error_sqr(tvi, ctoa)
                    if not self.is_wideband()
                    else vl.scaled_toa_error_sqr(tvi.toa, ctoa.toa_correction)
                )
                for (tvi, ctoa) in zip(self.toas, ctoas)
            ]
        )

    def model_dm(self, params: np.ndarray) -> np.ndarray:
        """Compute the model DM (dmu) for a given set of parameters."""
        params = vl.read_params(self.model, params)

        if not self.is_wideband():
            dms = np.zeros(len(self.toas))
            for ii, toa in enumerate(self.toas):
                ctoa = vl.TOACorrection()
                for component in self.model.components:
                    if vl.isa(component, vl.DispersionComponent):
                        dms[ii] += vl.value(
                            vl.dispersion_slope(component, toa, ctoa, params)
                        )
                    ctoa = vl.correct_toa(component, toa, ctoa, params)
        else:
            ctoas = [vl.correct_toa(self.model, tvi, params) for tvi in self.toas]
            dms = np.array([vl.value(ctoa.dm_correction.model_dm) for ctoa in ctoas])

        dmu_conversion_factor = 2.41e-16  # Hz / (DMconst * dmu)

        return dms * dmu_conversion_factor

    @classmethod
    def load_jlso(cls, filename: str):
        """Construct an `SPNTA` object from a JLSO file"""
        spnta = cls.__new__(cls)
        model, toas = vl.load_pulsar_data(filename)
        spnta.jlsofile = filename
        spnta._setup(model, toas)
        return spnta

    @classmethod
    def from_pint(
        cls, model: TimingModel, toas: TOAs, cheat_prior_scale=20, custom_priors={}
    ):
        """Construct an `SPNTA` object from PINT `TimingModel` and `TOAs` objects"""
        spnta = cls.__new__(cls)

        setup_log(level="WARNING")

        spnta.model_pint = model

        # custom_priors_dict is in the "raw" format. The numbers may be
        # in "normal" units and have to be converted into internal units.
        if isinstance(custom_priors, dict):
            custom_priors_dict = custom_priors
        elif isinstance(custom_priors, str):
            with open(custom_priors) as custom_priors_file:
                custom_priors_dict = json.load(custom_priors_file)
        else:
            custom_priors_dict = json.load(custom_priors)

        custom_priors = process_custom_priors(custom_priors_dict, model)

        model_v, toas_v = convert_model_and_toas(
            model,
            toas,
            cheat_prior_scale=cheat_prior_scale,
            custom_priors=custom_priors,
        )

        spnta._setup(model_v, toas_v)

        return spnta

    def update_pint_model(self, samples: np.ndarray) -> TimingModel:
        """Return an updataed PINT `TimingModel` based on posterior samples."""
        mp: TimingModel = deepcopy(self.model_pint)

        scaled_samples = self.rescale_samples(samples)

        for ii, pname in enumerate(self.param_names):
            if pname in mp.free_params:
                param_val = np.mean(scaled_samples[:, ii])
                param_err = np.std(scaled_samples[:, ii])
                mp[pname].value = param_val
                mp[pname].uncertainty_value = param_err

        return mp
