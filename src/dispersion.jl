export DispersionTaylor, dispersion_slope

"""Taylor series representation of the dispersion measure."""
struct DispersionTaylor <: DispersionComponent end

"""Compute the dispersion slope corresponding to a TOA."""
function dispersion_slope(::DispersionTaylor, ctoa::CorrectedTOA, params)
    t0 = params.DMEPOCH
    t = corrected_toa_value(ctoa)
    dms = params.DM
    dm = taylor_horner(t - t0, dms)
    return dm
end

const priors_DM = (
    LogUniform(1e+14, 1e+20),
    Uniform(-1e+16, 1e+16),
    Uniform(-1e+14, 1e+14),
    Uniform(-1e+14, 1e+14),
    Uniform(-1e+14, 1e+14),
)

"""
Prior distributions for DM ... DM4 based on psrcat values, but broadened to 
include one extra order of magnitude in either direction.
"""
lnprior(::DispersionTaylor, params::NamedTuple) = _lnprior(priors_DM, params, :DM)
