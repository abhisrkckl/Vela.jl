export ChromaticComponent, ChromaticTaylor, ChromaticPiecewise, chromatic_slope

"""
    ChromaticComponent

Abstract base type of all timing model components which provide a chromatic measure correction."""
abstract type ChromaticComponent <: DelayComponent end

"""Compute a chromatic delay."""
delay(component::ChromaticComponent, toa::TOA, toacorr::TOACorrection, params::NamedTuple) =
    chromatic_slope(component, toa, toacorr, params) *
    (frequency(1e6) / doppler_corrected_observing_frequency(toa, toacorr))^params.TNCHROMIDX

"""
    ChromaticTaylor
    
Taylor series representation of the chromatic measure.

Corresponds to `ChromaticCM` in `PINT`."""
struct ChromaticTaylor <: ChromaticComponent end

"""Compute the chromatic slope corresponding to a TOA using a Taylor series representation."""
function chromatic_slope(
    ::ChromaticTaylor,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)
    t0 = params.CMEPOCH
    t = corrected_toa_value(toa, toacorr, Float64)
    cms = params.CM
    cm = taylor_horner(t - t0, cms)
    return cm
end

"""Piecewise-constant representation of the chromatic measure with a 
constant chromatic index.

Corresponds to `ChromaticCMX` in `PINT`. Does not support overlapping 
ranges.
"""
struct ChromaticPiecewise <: ChromaticComponent
    cmx_mask::Vector{UInt}
end

function chromatic_slope(
    cmx::ChromaticPiecewise,
    toa::TOA,
    ::TOACorrection,
    params::NamedTuple,
)
    if is_tzr(toa)
        return zero(params.CMX_[1])
    end
    cmx_idx = cmx.cmx_mask[toa.index]
    return (cmx_idx == 0) ? zero(params.CMX_[1]) : params.CMX_[cmx_idx]
end
