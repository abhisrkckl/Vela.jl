from typing import IO, Optional
import emcee
import numpy as np
from scipy.optimize import dual_annealing
from numdifftools import Hessian
from pint.fitter import Fitter
from pint.models import TimingModel
from pint.toa import TOAs
from pint.residuals import Residuals, WidebandTOAResiduals
from pint.logging import log

from .spnta import SPNTA
from .vela import vl


class VelaFitter(Fitter):
    """An implementation of the PINT `Fitter` interface using `Vela` as a backend,
    applicable to both narrowband and wideband TOAs. This is a thin wrapper over the
    `SPNTA` class in `pyvela`.

    The fitting is done by sampling the posterior distribution using `emcee` or by maximizing
    the log-posterior using the `scipy.optimize.dual_annealing()` function. The former is more
    robust whereas the latter is faster. In the latter case, the uncertainties are estimated
    using the numerical Hessian at the maximum-posterior point. This may be unreliable when the
    posterior distribution is highly non-Gaussian. The fitting method can be chosen using the `mcmc`
    argument in the `fit_toas()` method.

    Note that unlike other `Fitter`s, the `model` and `toas` objects herein are immutable, and
    changing their contents will not make any difference in the fitting. If `toas` and `model`
    must be changed, please create a new instance of this class.

    Unlike the other fitters, `VelaFitter` needs the prior distributions of the free
    parameters to be specified, except for those with default priors. For other
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
            marginalize_gp_noise=True,
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
        """The PINT `TimingModel` object. Mutating this object has no effect on the fit.
        If this must be updated, create a new `VelaFitter` instance."""
        return self.spnta.model_pint

    @property
    def toas(self):
        """The PINT `TOAs` object. Mutating this object has no effect on the fit.
        If this must be updated, create a new `VelaFitter` instance."""
        return self.spnta.toas_pint

    def fit_toas_mcmc(
        self, nsteps: int = 6000, burnin: int = 2000, thin: int = 100, **kwargs
    ):
        """Fit the model to data. The parameter estimation is done by drawing samples from
        the posterior distribution using `emcee`.

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

        param_uncertainties = np.std(samples, axis=0)
        params_median = np.median(samples, axis=0)

        F0_ = np.longdouble(self.spnta.model.param_handler._default_params_tuple.F_.x)
        params_offset = np.array(
            [F0_ * int(par == "F0") for par in self.spnta.param_names]
        )

        return params_median + params_offset, param_uncertainties, samples

    def fit_toas_maxpost(
        self, maxiter: int = 1000, hessian_step: float = 1e-6, **kwargs
    ):
        """Fit the model to data. The best-fit parameters are estimated by maximizing the
        the log-posterior using `scipy.optimize.dual_annealing`. The parameter uncertainties
        are estimated using the numerical Hessian of the log-posterior at the maximum point. The
        uncertainties may not be accurate when the posterior distribution is highly non-Gaussian.

        Parameters
        ----------
        maxiter: int
            Number of dual annealing iterations
        hessian_step: int
            Step size used for computing the numerical Hessian
        """

        def _mlnpostq(cube: np.ndarray) -> float:
            """Compute the negative of the log-posterior given a point from the unit hypercube.
            The unit hypercube point is converted into a parameter space point using the
            `prior_transform`. This takes care of the"""
            return -self.spnta.lnpost(self.spnta.prior_transform(cube))

        cube0 = [0.5] * self.spnta.ndim
        bounds = [(0, 1)] * self.spnta.ndim
        result = dual_annealing(_mlnpostq, bounds, maxiter=maxiter, x0=np.array(cube0))

        params_raw = self.spnta.prior_transform(result.x)
        F0_ = np.longdouble(self.spnta.model.param_handler._default_params_tuple.F_.x)
        params_offset = np.array(
            [F0_ * int(par == "F0") for par in self.spnta.param_names]
        )
        params = self.spnta.rescale_samples(params_raw) + params_offset

        # If the a maximum-posterior value hits the bounds, adjust it slightly so that
        # we can compute the Hessian.
        dcube = (
            2
            * hessian_step
            * ((result.x == 0).astype(int) - (result.x == 1).astype(int))
        )

        cube1 = result.x + dcube
        H = Hessian(_mlnpostq, step=hessian_step)(cube1)
        assert np.all(np.isfinite(H)), "The Hessian contains non-finite elements!"

        epsilon_cube = np.sqrt(np.diag(np.linalg.pinv(H)))

        params1 = self.spnta.prior_transform(cube1)
        epsilon_params_raw = np.array(
            [
                eps_q / vl.pdf(prior.distribution, a)
                for a, eps_q, prior in zip(
                    params1, epsilon_cube, self.spnta.model.priors
                )
            ]
        )
        epsilon_params = self.spnta.rescale_samples(epsilon_params_raw)

        return params, epsilon_params, None

    def fit_toas(self, mcmc: bool = False, **kwargs):
        """Fit the model to data. The parameter estimation can be done by either sampling the
        posterior distribution using MCMC or by maximizing the log-posterior. The former is 
        more robust but the latter is faster. 

        The fitting methods are implemented in `fit_toas_mcmc()` and `fit_toas_maxpost()` methods.
        This function is just a wrapper around them.
        
        Parameters
        ----------
        mcmc: bool
            Whether to run MCMC or maximization for parameter estimation.
        **kwargs:
            Arguments forwarded to `fit_toas_mcmc()` or `fit_toas_maxpost()`.
            See their docstrings for more details.
        """
        params, uncertainties, self.samples = (
            self.fit_toas_mcmc(**kwargs) if mcmc else self.fit_toas_maxpost(**kwargs)
        )

        self.model.set_param_values(dict(zip(self.spnta.param_names, params)))
        self.model.set_param_uncertainties(
            dict(zip(self.spnta.param_names, uncertainties))
        )

        if not self.toas.wideband:
            self.resids.update()
        else:
            self.resids.toa.update()
            self.resids.dm.update()
