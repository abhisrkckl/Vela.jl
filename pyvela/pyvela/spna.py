from functools import cached_property
from typing import Iterable, List

import numpy as np

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
