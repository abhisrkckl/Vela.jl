from typing import IO, Optional
import emcee
import numpy as np
from pint.fitter import Fitter
from pint.models import TimingModel
from pint.toa import TOAs
from pint.residuals import Residuals, WidebandTOAResiduals
from pint.logging import log

from .spnta import SPNTA


class VelaFitter(Fitter):
    """An implementation of the PINT `Fitter` interface using `Vela` as a backend,
    applicable to both narrowband and wideband TOAs. This is a thin wrapper over the
    `SPNTA` class in `pyvela`. The fitting is done by sampling the posterior distribution
    using `emcee`.

    Note that unlike other `Fitter`s, the `toas` object herein is immutable, and changing its
    contents will not make any difference in the fitting.

    Unlike the other fitters, `VelaFitter` needs the prior distributions of the free
    parameters to be specified, except for those which have default priors. For other
    free parameters, either an `uncertainty_value` should be set from which a 'cheat'
    prior can be constructed, or a custom prior should be given through the `custom_priors`
    option. `custom_priors` can be a `dict` or a filename / IO object containing the prior
    in JSON format.

    The `prior` attribute in the PINT `Parameter` objects are *NOT* used, unlike `MCMCFitter`,
    since they are too slow.
    """

    def __init__(
        self,
        toas: TOAs,
        model: TimingModel,
        cheat_prior_scale: float = 100.0,
        custom_priors: dict | str | IO = {},
    ):
        log.info(
            "Constructing the `VelaFitter`... This will take some time due to JIT compilation."
        )

        self.spnta: SPNTA = SPNTA.from_pint(
            model,
            toas,
            cheat_prior_scale=cheat_prior_scale,
            custom_priors=custom_priors,
        )
        self.resids = (
            Residuals(self.toas, self.model)
            if not self.toas.wideband
            else WidebandTOAResiduals(self.toas, self.model)
        )

        self.samples: Optional[np.ndarray] = None

    @property
    def model(self):
        return self.spnta.model_pint

    @property
    def toas(self):
        return self.spnta.toas_pint

    def fit_toas_mcmc(
        self, nsteps: int = 6000, burnin: int = 2000, thin: int = 100, **kwargs
    ):
        """Fit the model to data.

        Parameters
        ----------
        nsteps: int
            Number of ensemble MCMC iterations
        burnin: int
            Number of samples to be discarded for burnin
        thin: int
            Thinning factor for MCMC samples
        """
        nwalkers = self.spnta.ndim * 5
        p0 = np.array(
            [
                self.spnta.prior_transform(cube)
                for cube in np.random.rand(nwalkers, self.spnta.ndim)
            ]
        )

        sampler = emcee.EnsembleSampler(
            nwalkers,
            self.spnta.ndim,
            self.spnta.lnpost_vectorized,
            vectorize=True,
            backend=emcee.backends.HDFBackend(f"__{self.model["PSR"].value}_chain.h5"),
        )

        sampler.run_mcmc(p0, nsteps, progress=True, progress_kwargs={"mininterval": 1})
        samples_raw = sampler.get_chain(flat=True, discard=burnin, thin=thin)
        samples = self.spnta.rescale_samples(samples_raw)

        self.samples = samples

        param_uncertainties = np.std(samples, axis=0)
        params_median = np.median(samples, axis=0)

        F0_ = np.longdouble(self.spnta.model.param_handler._default_params_tuple.F_.x)
        params_offset = np.array(
            [F0_ * int(par == "F0") for par in self.spnta.param_names]
        )

        self.model.set_param_values(
            dict(zip(self.spnta.param_names, params_median + params_offset))
        )
        self.model.set_param_uncertainties(
            dict(zip(self.spnta.param_names, param_uncertainties))
        )

        if not self.toas.wideband:
            self.resids.update()
        else:
            self.resids.toa.update()

    def fit_toas(self, mcmc:bool=False, **kwargs):
        if mcmc:
            self.fit_toas_mcmc(**kwargs)
        else:
            raise NotImplementedError
            # self.fit_toas_maxpost(**kwargs)
