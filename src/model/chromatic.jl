export ChromaticComponent, ChromaticTaylor, chromatic_slope

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
