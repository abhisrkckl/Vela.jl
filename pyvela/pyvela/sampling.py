import numpy as np

from .spnta import SPNTA


def get_start_samples(spnta: SPNTA, s: float, nwalkers: int) -> np.ndarray:
    """Get starting samples for the MCMC. nwalkers is the number of samples
    to be returned."""

    p0_ = np.array(
        [spnta.prior_transform(cube) for cube in np.random.rand(nwalkers, spnta.ndim)]
    )
    p0 = (
        ((1 - s) * spnta.maxpost_params + s * p0_)
        if np.isfinite(spnta.lnpost(spnta.default_params))
        else p0_
    )
    return p0
