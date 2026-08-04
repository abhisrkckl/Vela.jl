from functools import cached_property
from typing import Iterable, List, Tuple

import numpy as np
from scipy.linalg import cho_solve, cholesky, solve_triangular

from .spnta import SPNTA
from .vela import vl, jl


class SPNA:
    def __init__(self, spnta: SPNTA):
        self.spnta = spnta
        self.spna = vl.make_SPNA(
            spnta.model, spnta.toas, jl.Matrix(self.timing_model_designmatrix)
        )

    def lnpost(self, params):
        return vl.calc_lnpost(self.spna, params)

    def lnpost_vectorized(self, paramss):
        return vl.calc_lnpost_vectorized(self.spna, paramss)

    @cached_property
    def timing_model_designmatrix(self) -> np.ndarray:
        spnta = self.spnta
        M = np.empty((len(spnta.toas), spnta.ntmdim))
        F0 = spnta.model_pint_modified["F0"].quantity.to_value("Hz")
        dt = spnta.model_pint_modified.delay(spnta.toas_pint)
        for ii, (pname, sfac) in enumerate(
            zip(spnta.param_names[: spnta.ntmdim], spnta.scale_factors[: spnta.ntmdim])
        ):
            M[:, ii] = (
                spnta.model_pint_modified.d_phase_d_param(
                    spnta.toas_pint, dt, pname
                ).value
                / F0
                / sfac
            ).astype(float)
        return M

    def prior_transform(self, cube: Iterable[float]) -> Iterable[float]:
        """Compute the prior transform"""
        return vl.prior_transform(self.spna.priors, cube)

    def draw_from_prior(self, size: int = 1) -> np.ndarray:
        """Draw samples from prior."""
        return np.array(
            [self.prior_transform(np.random.rand(self.ndim)) for _ in range(size)]
        )

    @cached_property
    def param_names(self) -> Iterable[str]:
        """Free parameter names in the correct order. The names are same in both `Vela` and `PINT`,
        but the order may be different."""
        return np.array(list(vl.get_free_param_names(self.spna.param_handler)))

    @cached_property
    def ndim(self) -> int:
        """Number of free parameters."""
        return len(self.param_names)

    @cached_property
    def default_params(self) -> Iterable[str]:
        """Default parameter values taken from the par file."""
        return np.array(vl.read_param_values_to_vector(self.spna.param_handler))

    @cached_property
    def scale_factors(self) -> Iterable[float]:
        """Scale factors for converting free parameters from `PINT` units to `Vela` units."""
        return np.array(vl.get_scale_factors(self.spna.param_handler))

    @cached_property
    def marginalized_param_names(self) -> List[str]:
        """List of analytically marginalized parameters."""
        return list(vl.get_marginalized_param_names(self.spna))

    @cached_property
    def marginalized_default_params(self) -> np.ndarray:
        """Default values of analytically marginalized parameters."""
        pdict = dict(zip(self.spnta.param_names, self.spnta.default_params)) | dict(
            zip(
                self.spnta.marginalized_param_names,
                self.spnta.marginalized_default_params,
            )
        )
        return np.array([pdict[pname] for pname in self.marginalized_param_names])

    @cached_property
    def marginalized_param_scale_factors(self) -> np.ndarray:
        """Unit conversion factors for analytically marginalized parameters."""
        pdict = dict(zip(self.spnta.param_names, self.spnta.scale_factors)) | dict(
            zip(
                self.spnta.marginalized_param_names,
                self.spnta.marginalized_param_scale_factors,
            )
        )
        return np.array([pdict[pname] for pname in self.marginalized_param_names])

    def get_marginalized_param_offset_mean_and_covinvcf(
        self, params: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Returns the mean offsets and the inverse-covariance matrix Cholesky factor of the
        analytically marginalized parameters. Offsets are defined w.r.t. the default values.
        """
        params_ = vl.read_params(self.spna.param_handler, params)
        y, Ninvdiag = vl._calc_resids_and_Ninvdiag(self.spna, params_)
        M = np.array(self.spna.kernel.noise_basis)
        Phiinv = np.array(vl.calc_noise_weights_inv(self.spna.kernel, params_))
        Ninv_M = M * np.array(Ninvdiag)[:, None]
        MT_Ninv_y = y @ Ninv_M
        Sigmainv = np.diag(Phiinv) + M.T @ Ninv_M
        Sigmainv_cf = cholesky(Sigmainv, lower=False)
        ahat = cho_solve((Sigmainv_cf, False), MT_Ninv_y)
        return ahat, Sigmainv_cf

    def get_marginalized_param_offset_sample(
        self, params: np.ndarray
    ) -> Tuple[np.ndarray, float]:
        """Draw a sample of the analytically marginalized parameter vector given other parameters."""
        ahat, Sigmainv_cf = self.get_marginalized_param_offset_mean_and_covinvcf(params)
        z = np.random.randn(len(ahat))
        da = solve_triangular(Sigmainv_cf, z, lower=False)
        a = ahat + da
        lnp = (
            -0.5 * (z @ z)
            + np.sum(np.log(np.diag(Sigmainv_cf)))
            - 0.5 * len(ahat) * np.log(2 * np.pi)
        )

        return a, lnp

    def get_marginalized_param_sample(
        self, params: np.ndarray
    ) -> Tuple[np.ndarray, float]:
        """Draw a sample of the analytically marginalized parameter values given other parameters."""
        da, lnp = self.get_marginalized_param_offset_sample(params)
        return self.marginalized_default_params + da, lnp

    @cached_property
    def ntmdim(self):
        return self.spnta.ntmdim

    def get_spnta_sample(self, params: np.ndarray):
        mpar, lnp_cond = self.get_marginalized_param_sample(params)
        lnp_spna = self.lnpost(params)
        return np.hstack((mpar[: self.ntmdim], params)), lnp_spna + lnp_cond
