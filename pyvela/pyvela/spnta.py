import datetime
from functools import cached_property
import getpass
import json
from copy import deepcopy
import os
import platform
import sys
from typing import IO, Dict, Iterable, List, Optional
import warnings

import astropy
import emcee
import numpy as np
import pint
import pyvela
from scipy.linalg import cho_factor, cho_solve
from pint.binaryconvert import convert_binary
from pint.logging import setup as setup_log
from pint.models import TimingModel, get_model, get_model_and_toas
from pint.toa import TOAs

from .ecorr import ecorr_sort
from .model import fix_params, pint_model_to_vela
from .parameters import get_unit_conversion_factor
from .priors import process_custom_priors
from .toas import day_to_s, pint_toas_to_vela
from .vela import jl, vl


def convert_model_and_toas(
    model: TimingModel,
    toas: TOAs,
    noise_params: List[str],
    marginalize_gp_noise: bool,
    cheat_prior_scale: float = 100.0,
    custom_priors: dict = {},
):
    """Read a pair of par & tim files and create a `Vela.TimingModel` object and a
    Julia `Vector` of `TOA`s."""

    fix_params(model, toas)

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
        noise_params,
        marginalize_gp_noise,
        ecorr_toa_ranges=ecorr_toa_ranges,
        ecorr_indices=ecorr_indices,
    )
    toas_v = pint_toas_to_vela(toas, float(model["PEPOCH"].value))

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
    param_names: array-like[str]
        Free parameter names
    param_units: array-like[str]
        Free parameter units
    param_labels: array-like[str]
        Free parameter labels (for plotting)
    scale_factors: ndarray
        Scale factors for converting between Vela's internal units and the commonly used units.
    maxlike_params: ndarray
        Free parameter values taken from the par file in internal units
    ndim: int
        Number of free parameters
    """

    def __init__(
        self,
        parfile: str,
        timfile: str,
        marginalize_gp_noise: bool = True,
        cheat_prior_scale: float = 100.0,
        custom_priors: str | IO | dict = {},
        check: bool = True,
        pint_kwargs: dict = {},
    ):
        self.parfile = parfile
        self.timfile = timfile
        self.jlsofile: Optional[str] = None
        self.custom_prior_file: Optional[str] = None

        self.cheat_prior_scale: Optional[float] = cheat_prior_scale

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
        self.model_pint = deepcopy(model_pint)
        self.model_pint_modified = model_pint
        self.toas_pint = toas_pint

        # custom_priors_dict is in the "raw" format. The numbers may be
        # in "normal" units and have to be converted into internal units.
        if isinstance(custom_priors, dict):
            self.custom_priors_dict = custom_priors
        elif isinstance(custom_priors, str):
            self.custom_prior_file = custom_priors
            with open(custom_priors) as custom_priors_file:
                self.custom_priors_dict = json.load(custom_priors_file)
        else:
            self.custom_priors_dict = json.load(custom_priors)

        custom_priors = process_custom_priors(self.custom_priors_dict, model_pint)

        # Use the original PINT TimingModel object.
        noise_params = self.model_pint.get_params_of_component_type("NoiseComponent")

        setup_log(level="WARNING")
        model, toas = convert_model_and_toas(
            model_pint,
            toas_pint,
            noise_params,
            marginalize_gp_noise,
            cheat_prior_scale=cheat_prior_scale,
            custom_priors=custom_priors,
        )

        self.pulsar = vl.Pulsar(model, toas)

        if check:
            self._check()

    def _check(self):
        """Check if the computations work with the default values."""
        cube = np.random.rand(self.ndim)
        sample = self.prior_transform(cube)
        lnpr = self.lnprior(sample)
        lnl = self.lnlike(sample)
        lnp = self.lnpost(sample)
        lnpv = self.lnpost_vectorized(np.array([sample]))
        assert all(np.isfinite([lnpr, lnl, lnp])), (
            "The log-prior, log-likelihood, or log-posterior is non-finite at the default parameter values. "
            "Please make sure that (1) default values are within the prior range, and (b) the default values "
            "provide a phase-connected solution. If nothing is wrong, this may be a bug. Please report this."
        )
        assert np.isclose(lnp, lnpv), (
            "The non-vectorized and vectorized versions of the log-posterior gives different results. "
            "This is most likely a bug. Please report this."
        )

    def lnlike(self, params: Iterable[float]) -> float:
        """Compute the log-likelihood function"""
        return vl.calc_lnlike(self.pulsar, params)

    def lnprior(self, params: Iterable[float]) -> float:
        """Compute the log-prior distribution"""
        return vl.calc_lnprior(self.pulsar, params)

    def prior_transform(self, cube: Iterable[float]) -> Iterable[float]:
        """Compute the prior transform"""
        return vl.prior_transform(self.pulsar, cube)

    def lnpost(self, params: Iterable[float]) -> float:
        """Compute the log-posterior distribution"""
        return vl.calc_lnpost(self.pulsar, params)

    def lnpost_vectorized(self, paramss: np.ndarray) -> Iterable[float]:
        """Compute the log-posterior distribution over a collection
        of points in the parameter space"""
        return vl.calc_lnpost_vectorized(self.pulsar, paramss)

    @property
    def model(self):
        """The `Vela.TimingModel` object."""
        return self.pulsar.model

    @property
    def toas(self):
        """The `Vector{Vela.TOA}` or `Vector{Vela.WidebandTOA}` object."""
        return self.pulsar.toas

    @cached_property
    def param_names(self) -> Iterable[str]:
        """Free parameter names in the correct order. The names are same in both `Vela` and `PINT`,
        but the order may be different."""
        return np.array(list(vl.get_free_param_names(self.pulsar.model)))

    @cached_property
    def param_labels(self) -> Iterable[str]:
        """Free parameter labels containing parameter names and units."""
        return np.array(list(vl.get_free_param_labels(self.pulsar.model)))

    @cached_property
    def param_units(self) -> Iterable[str]:
        """String representations of `PINT` units of free parameters. Tfhese strings are supported by
        `astropy.units.`"""
        return np.array(list(vl.get_free_param_units(self.pulsar.model)))

    @cached_property
    def param_prefixes(self) -> Iterable[str]:
        """Free parameter prefixes. For non-mask/prefix parameters the prefix is the same as
        the parameter name."""
        return np.array(list(vl.get_free_param_prefixes(self.pulsar.model)))

    @cached_property
    def scale_factors(self) -> Iterable[float]:
        """Scale factors for converting free parameters from `PINT` units to `Vela` units."""
        return np.array(vl.get_scale_factors(self.pulsar.model))

    @cached_property
    def default_params(self) -> Iterable[str]:
        """Default parameter values taken from the par file."""
        return np.array(vl.read_param_values_to_vector(self.pulsar.model))

    @cached_property
    def ndim(self) -> int:
        """Number of free parameters."""
        return len(self.param_names)

    @cached_property
    def ntmdim(self) -> int:
        """Number of free timing model parameters (does not include noise parameters)."""
        return vl.get_num_timing_params(self.pulsar.model)

    @cached_property
    def has_marginalized_gp_noise(self) -> bool:
        """Whether the model contains marginalized correlated Gaussian noise processes."""
        return vl.isa(self.model.kernel, vl.WoodburyKernel)

    @cached_property
    def has_ecorr_noise(self) -> bool:
        return vl.isa(self.model.kernel, vl.EcorrKernel) or (
            self.has_marginalized_gp_noise
            and vl.isa(self.model.kernel.inner_kernel, vl.EcorrKernel)
        )

    def get_marginalized_gp_noise_realization(self, params: np.ndarray) -> np.ndarray:
        """Get a realization of the marginalized GP noise given a set of parameters.
        The length of `params` should be the same as the number of free parameters."""
        assert self.has_marginalized_gp_noise
        params_ = vl.read_params(self.model, params)
        y, Ndiag = vl._calc_resids_and_Ndiag(self.model, self.toas, params_)
        M = np.array(self.model.kernel.noise_basis)
        Phiinv = np.array(vl.calc_noise_weights_inv(self.model.kernel, params_))
        Ninv_M = M / np.array(Ndiag)[:, None]
        MT_Ninv_y = y @ Ninv_M
        Sigmainv = np.diag(Phiinv) + M.T @ Ninv_M
        Sigmainv_cf = cho_factor(Sigmainv)
        ahat = cho_solve(Sigmainv_cf, MT_Ninv_y)
        return M @ ahat

    def rescale_samples(self, samples: np.ndarray) -> np.ndarray:
        """Rescale the samples from Vela's internal units to common units"""
        return samples / self.scale_factors

    def save_jlso(self, filename: str) -> None:
        """Write the model and TOAs as a JLSO file"""
        vl.save_pulsar_data(filename, self.pulsar.model, self.pulsar.toas)

    @cached_property
    def wideband(self) -> bool:
        """Whether the TOAs are wideband."""
        return jl.isa(self.pulsar.toas[0], vl.WidebandTOA)

    @cached_property
    def mjds(self) -> np.ndarray:
        """Get the MJDs of each TOA."""
        if self.wideband:
            return np.array(
                [
                    jl.Float64(vl.value(wtoa.toa.value)) / day_to_s
                    for wtoa in self.pulsar.toas
                ]
            )
        else:
            return np.array(
                [jl.Float64(vl.value(toa.value)) / day_to_s for toa in self.pulsar.toas]
            )

    def time_residuals(self, params: np.ndarray) -> np.ndarray:
        """Get the timing residuals (s) for a given set of parameters."""
        params = vl.read_params(self.pulsar.model, params)
        return np.array(
            [
                vl.value(wr[0])
                for wr in vl.form_residuals(self.pulsar.model, self.pulsar.toas, params)
            ]
            if self.wideband
            else list(
                map(
                    vl.value,
                    vl.form_residuals(self.pulsar.model, self.pulsar.toas, params),
                )
            )
        )

    def whitened_time_residuals(self, params: np.ndarray) -> np.ndarray:
        """Get whitened time residuals using the given set of parameters. This is done by
        subtracting the marginalized GP noise realizations from the time residuals."""
        return (
            self.time_residuals(params)
            - self.get_marginalized_gp_noise_realization(params)[: len(self.toas)]
            if self.has_marginalized_gp_noise
            else self.time_residuals(params)
        )

    def dm_residuals(self, params: np.ndarray) -> np.ndarray:
        """Get the DM residuals (s) for a given set of parameters (wideband only)."""
        assert self.wideband, "This method is only defined for wideband datasets."

        params = vl.read_params(self.pulsar.model, params)

        return np.array(
            [
                vl.value(wr[1])
                for wr in vl.form_residuals(self.pulsar.model, self.pulsar.toas, params)
            ]
        )

    def whitened_dm_residuals(self, params: np.ndarray) -> np.ndarray:
        """Get whitened DM residuals using the given set of parameters. This is done by
        subtracting the marginalized GP noise realizations from the DM residuals."""
        return (
            self.dm_residuals(params)
            - self.get_marginalized_gp_noise_realization(params)[len(self.toas) :]
            if self.has_marginalized_gp_noise
            else self.dm_residuals(params)
        )

    def scaled_toa_unceritainties(self, params: np.ndarray) -> np.ndarray:
        """Get the scaled TOA uncertainties (s) for a given set of parameters."""
        params = vl.read_params(self.pulsar.model, params)

        ctoas = [
            vl.correct_toa(self.pulsar.model, tvi, params) for tvi in self.pulsar.toas
        ]

        return np.sqrt(
            [
                vl.value(
                    vl.scaled_toa_error_sqr(tvi.toa, ctoa.toa_correction)
                    if self.wideband
                    else vl.scaled_toa_error_sqr(tvi, ctoa)
                )
                for (tvi, ctoa) in zip(self.pulsar.toas, ctoas)
            ]
        )

    def scaled_dm_unceritainties(self, params: np.ndarray) -> np.ndarray:
        """Get the scaled DM uncertainties (s) for a given set of parameters (wideband only)."""
        assert self.wideband, "This method is only defined for wideband datasets."

        params = vl.read_params(self.pulsar.model, params)

        ctoas = [
            vl.correct_toa(self.pulsar.model, tvi, params) for tvi in self.pulsar.toas
        ]

        return np.sqrt(
            [
                vl.value(vl.scaled_dm_error_sqr(tvi.dminfo, ctoa.dm_correction))
                for (tvi, ctoa) in zip(self.pulsar.toas, ctoas)
            ]
        )

    def model_dm(self, params: np.ndarray) -> np.ndarray:
        """Compute the model DM (dmu) for a given set of parameters."""
        params = vl.read_params(self.pulsar.model, params)

        if not self.wideband:
            dms = np.zeros(len(self.pulsar.toas))
            for ii, toa in enumerate(self.pulsar.toas):
                ctoa = vl.TOACorrection()
                for component in self.pulsar.model.components:
                    if vl.isa(component, vl.DispersionComponent):
                        dms[ii] += vl.value(
                            vl.dispersion_slope(component, toa, ctoa, params)
                        )
                    ctoa = vl.correct_toa(component, toa, ctoa, params)
        else:
            ctoas = [
                vl.correct_toa(self.pulsar.model, tvi, params)
                for tvi in self.pulsar.toas
            ]
            dms = np.array([vl.value(ctoa.dm_correction.model_dm) for ctoa in ctoas])

        dmu_conversion_factor = 2.41e-16  # Hz / (DMconst * dmu)

        return dms * dmu_conversion_factor

    @classmethod
    def load_jlso(cls, jlsoname: str, parfile: str) -> "SPNTA":
        """Construct an `SPNTA` object from a JLSO file"""
        spnta = cls.__new__(cls)
        model, toas = vl.load_pulsar_data(jlsoname)

        spnta.jlsofile = jlsoname
        spnta.parfile = parfile
        spnta.timfile = None
        spnta.custom_prior_file = None
        spnta.cheat_prior_scale = None

        spnta.pulsar = vl.Pulsar(model, toas)
        spnta.model_pint = get_model(parfile)
        spnta.model_pint_modified = None
        spnta.toas_pint = None
        spnta._check()
        return spnta

    @classmethod
    def from_pint(
        cls,
        model: TimingModel,
        toas: TOAs,
        marginalize_gp_noise: bool = True,
        cheat_prior_scale: float = 100.0,
        custom_priors: dict | str | IO = {},
    ) -> "SPNTA":
        """Construct an `SPNTA` object from PINT `TimingModel` and `TOAs` objects"""
        spnta = cls.__new__(cls)

        setup_log(level="WARNING")

        if not toas.planets:
            warnings.warn("Computing planetary ephemerides...")
            toas.compute_posvels(ephem=model["EPHEM"].value, planets=True)

        spnta.model_pint = deepcopy(model)
        spnta.model_pint_modified = model
        spnta.toas_pint = toas

        spnta.parfile = model.name
        spnta.timfile = toas.filename
        spnta.custom_prior_file = None
        spnta.jlsofile = None

        spnta.cheat_prior_scale = cheat_prior_scale

        # custom_priors_dict is in the "raw" format. The numbers may be
        # in "normal" units and have to be converted into internal units.
        if isinstance(custom_priors, dict):
            custom_priors_dict = custom_priors
        elif isinstance(custom_priors, str):
            spnta.custom_prior_file = custom_priors
            with open(custom_priors) as custom_priors_file:
                custom_priors_dict = json.load(custom_priors_file)
        else:
            custom_priors_dict = json.load(custom_priors)

        custom_priors = process_custom_priors(custom_priors_dict, model)

        # Use the original PINT TimingModel object.
        noise_params = spnta.model_pint.get_params_of_component_type("NoiseComponent")

        model_v, toas_v = convert_model_and_toas(
            model,
            toas,
            noise_params,
            marginalize_gp_noise,
            cheat_prior_scale=cheat_prior_scale,
            custom_priors=custom_priors,
        )

        spnta.pulsar = vl.Pulsar(model_v, toas_v)
        spnta._check()

        return spnta

    def update_pint_model(self, samples: np.ndarray) -> TimingModel:
        """Return an updated PINT `TimingModel` based on posterior samples."""
        mp: TimingModel = deepcopy(self.model_pint_modified)

        scaled_samples = self.rescale_samples(samples)

        for ii, pname in enumerate(self.param_names):
            if pname in mp.free_params:
                param_val = np.mean(scaled_samples[:, ii])
                param_err = np.std(scaled_samples[:, ii])
                mp[pname].value = param_val
                mp[pname].uncertainty_value = param_err

        return mp

    def full_prior_dict(self, scaled: bool = True):
        """Returns a dictionary containing prior information for all parameters."""
        result = {}
        for prior, pname, punit, scale_factor in zip(
            self.model.priors, self.param_names, self.param_units, self.scale_factors
        ):
            ptype = str(prior.source_type)
            dname = str(vl.distr_name(prior))
            dtype = (
                getattr(jl.Distributions, dname)
                if hasattr(jl.Distributions, dname)
                else getattr(vl, dname)
            )
            dargs = (
                vl.unscale_prior_args(
                    dtype,
                    vl.distr_args(prior),
                    scale_factor,
                )
                if hasattr(jl.Distributions, dname)
                else vl.distr_args(prior)
            )
            prior_dict = {
                pname: {
                    "distribution": dname,
                    "args": dargs,
                    "type": ptype,
                    "unit": punit,
                }
            }
            if jl.isa(prior.distribution, vl.Truncated):
                if prior.distribution.upper is not None:
                    prior_dict["upper"] = prior.distribution.upper / scale_factor
                if prior.distribution.lower is not None:
                    prior_dict["lower"] = prior.distribution.lower / scale_factor

            result.update(prior_dict)

        return result

    def info_dict(self, sampler_info: Dict = {}, truth_par_file: Optional[str] = None):
        """Returns a dictionary containing information about the machine, environment, sampler,
        input, etc."""
        info_dict = {
            "input": {
                "par_file": (
                    os.path.basename(self.parfile) if self.parfile is not None else None
                ),
                "tim_file": (
                    os.path.basename(self.timfile) if self.timfile is not None else None
                ),
                "jlso_file": (
                    os.path.basename(self.jlsofile)
                    if self.jlsofile is not None
                    else None
                ),
                "custom_prior_file": (
                    os.path.basename(self.custom_prior_file)
                    if self.custom_prior_file is not None
                    else None
                ),
                "cheat_prior_scale": self.cheat_prior_scale,
                "truth_par_file": (
                    os.path.basename(truth_par_file)
                    if truth_par_file is not None
                    else None
                ),
            },
            "sampler": sampler_info,
            "env": {
                "launch_time": datetime.datetime.now().isoformat(),
                "user": getpass.getuser(),
                "host": platform.node(),
                "os": platform.platform(),
                "julia_threads": vl.nthreads(),
                "python": sys.version,
                "julia": str(vl.VERSION),
                "pyvela": pyvela.__version__,
                "pint": pint.__version__,
                "emcee": emcee.__version__,
                "numpy": np.__version__,
                "astropy": astropy.__version__,
            },
        }

        return info_dict

    def save_new_parfile(
        self, params: np.ndarray, param_uncertainties: np.ndarray, filename: str
    ):
        """Save a new par file given a set of parameters and uncertainties."""
        param_vals = self.rescale_samples(params)
        param_errs = self.rescale_samples(param_uncertainties)

        model1 = (
            deepcopy(self.model_pint_modified)
            if self.model_pint_modified is not None
            else self.model_pint
        )
        for pname, pval, perr in zip(self.param_names, param_vals, param_errs):
            if pname in model1:
                model1[pname].value = (
                    pval
                    if pname != "F0"
                    else (
                        np.longdouble(
                            self.model.param_handler._default_params_tuple.F_.x
                        )
                        + pval
                    )
                )
                model1[pname].uncertainty_value = perr
            else:
                warnings.warn(
                    f"Parameter {pname} not found in the PINT TimingModel!"
                )  # pragma: no cover

        model1.write_parfile(filename)

    def save_resids(self, params: np.ndarray, filename: str) -> None:
        """Save the residuals and scaled uncertainties into a text file
        given a set of parameters."""
        wb = self.wideband

        ntoas = len(self.toas)
        mjds = self.mjds
        tres = self.time_residuals(params)
        tres_w = self.whitened_time_residuals(params)
        terr = self.scaled_toa_unceritainties(params)

        res_arr = np.zeros((ntoas, 3 * (1 + int(wb)) + 1))
        res_arr[:, 0] = mjds
        res_arr[:, 1] = tres
        res_arr[:, 2] = tres_w
        res_arr[:, 3] = terr

        if wb:
            dres = self.dm_residuals(params)
            dres_w = self.whitened_dm_residuals(params)
            derr = self.scaled_dm_unceritainties(params)

            res_arr[:, 4] = dres
            res_arr[:, 5] = dres_w
            res_arr[:, 6] = derr

        np.savetxt(filename, res_arr)

    def save_results(
        self,
        outdir: str,
        samples_raw: np.ndarray,
        sampler_info: dict = {},
        truth_par_file: Optional[str] = None,
    ) -> None:
        """Given the posterior samples, save the results into an output directory.
        `pyvela` script uses this function to save the results.

        `outdir` is the output directory to which the results will be saved.

        `samples_raw` is the burned-in and thinned MCMC chain obtained from the sampler.
        If these samples have associated importance weights (e.g., from nested sampling),
        please resample them before passing into this method such that each sample has
        equal weight.

        `sampler_info` is a dictionary containing sampler configuration information. It will
        be saved as-is into the summary file.

        `truth_par_file` is the original par file that was used to simulate a dataset. This is
        only applicable for simulated datasets.

        The following files are saved.

            1. `samples_raw.npy` - Samples in Vela's internal units (numpy format)
            2. `samples.npy` - Samples in 'normal' units (numpy format)
            3. `params_median.txt` - Parameter median values
            4. `params_std.txt` - Parameter standard deviations
            5. `<PSRNAME>.median.par` - par file containing median parameter values and standard deviations (PINT format)
            6. `residuals.txt` - Residuals (time for narrowband, time & DM for wideband) computed using median parameter values
            7. `param_default_values.txt` - Default parameter values taken from the input par file
            8. `param_names.txt` - Parameter names (PINT format)
            9. `param_prefixes.txt` - Parameter prefixes (PINT format)
            10. `param_units.txt` - Parameter unit strings (astropy format)
            11. `param_scale_factors.txt` - Scale factors that convert parameter values from Vela's internal units to 'normal' units
            12. `param_true_values.txt` - Parameter values used for simulation, taken from the 'truth' par file
            13. `param_autocorr.txt` - Parameter autocorrelation lengths (from the thinned chains)
            14. `prior_info.json` - Prior distributions on all free parameters (JSON format)
            15. `prior_evals.npy` - Prior distributions evaluated in the posterior range for plotting (numpy format)
            16. `summary.json` - Information about the machine, environment, sampler, and input (JSON format)
        """
        samples = self.rescale_samples(samples_raw)

        with open(f"{outdir}/samples_raw.npy", "wb") as f:
            np.save(f, samples_raw)
        with open(f"{outdir}/samples.npy", "wb") as f:
            np.save(f, samples)

        param_uncertainties = np.std(samples_raw, axis=0)
        params_median = np.median(samples_raw, axis=0)
        np.savetxt(f"{outdir}/params_median.txt", params_median)
        np.savetxt(f"{outdir}/params_std.txt", param_uncertainties)
        self.save_new_parfile(
            params_median,
            param_uncertainties,
            f"{outdir}/{self.model.pulsar_name}.median.par",
        )

        self.save_resids(params_median, f"{outdir}/residuals.txt")
        param_autocorr = emcee.autocorr.integrated_time(
            samples_raw, quiet=True, has_walkers=False
        )

        np.savetxt(f"{outdir}/param_default_values.txt", self.default_params)
        np.savetxt(f"{outdir}/param_names.txt", self.param_names, fmt="%s")
        np.savetxt(f"{outdir}/param_prefixes.txt", self.param_prefixes, fmt="%s")
        np.savetxt(f"{outdir}/param_units.txt", self.param_units, fmt="%s")
        np.savetxt(f"{outdir}/param_scale_factors.txt", self.scale_factors)
        np.savetxt(f"{outdir}/param_autocorr.txt", param_autocorr)

        if truth_par_file is not None:
            np.savetxt(
                f"{outdir}/param_true_values.txt", get_true_values(self, truth_par_file)
            )

        with open(f"{outdir}/prior_info.json", "w") as prior_info_file:
            json.dump(self.full_prior_dict(), prior_info_file, indent=4)

        self._save_prior_evals(samples, f"{outdir}/prior_evals.npy")

        summary_info = self.info_dict(sampler_info, truth_par_file)
        with open(f"{outdir}/summary.json", "w") as summary_file:
            json.dump(summary_info, summary_file, indent=4)

    def _single_param_prior(self, param_idx: int, value: float):
        prior = self.model.priors[param_idx]
        scale_factor = self.scale_factors[param_idx]
        return vl.pdf(prior.distribution, value * scale_factor)

    def _save_prior_evals(self, samples: np.ndarray, filename: str):
        """Save the prior PDF evaluated at uniformly spaced points within the
        posterior range."""
        nn = 1000
        result = np.empty((nn, 2 * self.ndim))

        for ii in range(self.ndim):
            xs = np.linspace(np.min(samples[:, ii]), np.max(samples[:, ii]), nn)
            ys = np.array([self._single_param_prior(ii, x) for x in xs])
            result[:, 2 * ii] = xs
            result[:, 2 * ii + 1] = ys

        np.save(filename, result)


def get_true_values(spnta: SPNTA, truth_par_file: str):
    """Read free parameter values from the "truth" par file containing
    original parameter values used for simulating a dataset."""
    true_model = get_model(truth_par_file)
    return (
        np.array(
            [
                (
                    (true_model[par].value if par != "F0" else 0.0)
                    if par in true_model
                    else np.nan
                )
                for par in spnta.param_names
            ]
        )
        * spnta.scale_factors
    )
