from typing import IO
import emcee
import numpy as np
from pint.fitter import Fitter
from pint.models import TimingModel
from pint.toa import TOAs
from pint.residuals import Residuals, WidebandTOAResiduals

from .spnta import SPNTA


class VelaFitter(Fitter):
    def __init__(
        self,
        toas: TOAs,
        model: TimingModel,
        cheat_prior_scale: float = 100.0,
        custom_priors: dict | str | IO = {},
    ):
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

    @property
    def model(self):
        return self.spnta.model_pint

    @property
    def toas(self):
        return self.spnta.toas_pint

    def fit_toas(
        self, nsteps: int = 6000, burnin: int = 1500, thin: int = 100, **kwargs
    ):
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
            backend=emcee.backends.HDFBackend("__chain.h5"),
        )

        sampler.run_mcmc(p0, nsteps, progress=True, progress_kwargs={"mininterval": 1})
        samples_raw = sampler.get_chain(flat=True, discard=burnin, thin=thin)
        samples = self.spnta.rescale_samples(samples_raw)

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
        self.resids.update()
