"""Abstrct base type of all timing model components which provide a chromatic measure correction."""
abstract type ChromaticComponent <: DelayComponent end

"""Compute a chromatic delay."""
delay(component::ChromaticComponent, ctoa::CorrectedTOA, params) =
    chromatic_slope(component, ctoa, params) *
    (frequency(1e6) / doppler_corrected_observing_frequency(ctoa))^params.TNCHROMIDX

"""Taylor series representation of the chromatic measure.

Corresponds to `DispersionDM` in `PINT`."""
struct ChromaticTaylor <: ChromaticComponent end

"""Compute the chromatic slope corresponding to a TOA using a Taylor series representation."""
function chromatic_slope(::ChromaticTaylor, ctoa::CorrectedTOA, params)
    t0 = params.CMEPOCH
    t = corrected_toa_value(ctoa)
    cms = params.CM
    cm = taylor_horner(t - t0, cms)
    return cm
end
