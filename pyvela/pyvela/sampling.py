import numpy as np

from .spnta import SPNTA


def get_start_samples(spnta: SPNTA, s: float, nwalkers: int) -> np.ndarray:
    """Get starting samples for the MCMC. nwalkers is the number of samples
    to be returned."""

    assert np.isfinite(
        spnta.lnpost(spnta.default_params)
    ), "The default values are outside the prior range. Please check your input par file."

    assert 0 < s <= 1

    p0_ = np.array(
        [spnta.prior_transform(cube) for cube in np.random.rand(nwalkers, spnta.ndim)]
    )

    p0 = np.empty_like(p0_)

    for ii, (pname, defval, sfac) in enumerate(
        zip(spnta.param_names, spnta.default_params, spnta.scale_factors)
    ):
        param = spnta.model_pint[pname] if pname in spnta.model_pint else None
        if (
            param is not None
            and param.uncertainty is not None
            and param.uncertainty_value > 0
        ):
            err = spnta.model_pint[pname].uncertainty_value * sfac
            p0[:, ii] = defval + err * 10 * s * np.random.randn(nwalkers)
        else:
            p0[:, ii] = (1 - s) * defval + s * p0_[:, ii]

    for w in range(nwalkers):
        lnp = spnta.lnprior(p0[w, :])
        if not np.isfinite(lnp):
            p0[w, :] = (1 - s) * spnta.default_params + s * p0_[w, :]

    return p0
