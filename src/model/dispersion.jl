export DispersionTaylor, dispersion_slope

"""Taylor series representation of the dispersion measure.

Corresponds to `DispersionDM` in `PINT`.

Reference:
    [Backer & Hellings 1986](http://doi.org/10.1146/annurev.aa.24.090186.002541)
"""
struct DispersionTaylor <: DispersionComponent end

function dispersion_slope(
    ::DispersionTaylor,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)::GQ
    t0 = params.DMEPOCH
    t = corrected_toa_value(toa, toacorr, Float64)
    dms = params.DM
    dm = taylor_horner(t - t0, dms)
    return dm
end
