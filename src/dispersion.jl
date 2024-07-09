export DispersionTaylor, dispersion_slope

"""Taylor series representation of the dispersion measure."""
struct DispersionTaylor <: DispersionComponent
    DMEPOCH::GQ{Float64}
end

"""Compute the dispersion slope corresponding to a TOA."""
function dispersion_slope(dmt::DispersionTaylor, ctoa::CorrectedTOA, params)
    t0 = dmt.DMEPOCH
    t = corrected_toa_value(ctoa)
    dms = params.DM
    dm = taylor_horner(t - t0, dms)
    return dm
end
