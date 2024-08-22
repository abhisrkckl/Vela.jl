export DispersionTaylor, dispersion_slope

"""Taylor series representation of the dispersion measure.

Corresponds to `DispersionDM` in `PINT`."""
struct DispersionTaylor <: DispersionComponent end

"""Compute the dispersion slope corresponding to a TOA using a Taylor series representation."""
function dispersion_slope(::DispersionTaylor, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    t0 = params.DMEPOCH
    t = corrected_toa_value(ctoa)
    dms = params.DM
    dm = taylor_horner(t - t0, dms)
    return dm
end
