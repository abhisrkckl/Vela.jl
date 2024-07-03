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

"""
Prior distributions for DM ... DM4 based on psrcat values, but broadened to 
include one extra order of magnitude in either direction.
"""
function lnprior(::DispersionTaylor, params::NamedTuple)
    priors_DM = (
        LogUniform(1e+14, 1e+20),
        Uniform(-1e+16, 1e+16),
        Uniform(-1e+14, 1e+14),
        Uniform(-1e+14, 1e+14),
        Uniform(-1e+14, 1e+14),
    )

    vals = map(value, params.DM)
    dists = priors_DM[1:length(params.DM)]

    return sum(map(logpdf, dists, vals))
end
