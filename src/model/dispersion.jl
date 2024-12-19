export DispersionTaylor, DispersionPiecewise, dispersion_slope

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

"""Piecewise-constant representation of the dispersion measure.

Corresponds to `DispersionDMX` in `PINT`.

Reference:
    [Arzoumanian+ 2015](http://doi.org/10.1088/0004-637X/813/1/65)
"""
struct DispersionPiecewise <: DispersionComponent
    dmx_mask::Vector{UInt}
end

function dispersion_slope(
    dmx::DispersionPiecewise,
    toa::TOA,
    ::TOACorrection,
    params::NamedTuple,
)
    dmx_idx = dmx.dmx_mask[toa.index]
    return (dmx_idx == 0) ? zero(params.DMX_[1]) : params.DMX_[dmx_idx]
end