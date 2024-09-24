export Component, PhaseComponent, DelayComponent, WhiteNoiseComponent, correct_toa

"""Abstract base type of all timing & noise model components which affect one TOA at a time."""
abstract type Component end

"""Abstract base type of all timing model components which contribute a phase correction to a TOA."""
abstract type PhaseComponent <: Component end

"""Correct the phase of a `TOA`."""
correct_toa(
    component::PhaseComponent,
    toa::TOABase,
    toacorr::TOACorrection,
    params::NamedTuple,
) = correct_toa_phase(toacorr; phase = phase(component, toa, toacorr, params))

"""Abstract base type of all timing model components which contribute a time delay correction to a TOA."""
abstract type DelayComponent <: Component end

"""Correct the value of a `TOA` using a delay."""
correct_toa(
    component::DelayComponent,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
) = correct_toa_delay(toacorr; delay = delay(component, toa, toacorr, params))

"""Abstrct base type of all timing model components which provide a dispersion measure correction."""
abstract type DispersionComponent <: DelayComponent end

"""Compute a dispersion delay."""
function delay(
    component::DispersionComponent,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)
    dm = dispersion_slope(component, toa, toacorr, params)
    ν = doppler_corrected_observing_frequency(toa, toacorr)
    return dm / (ν * ν)
end

"""Abstract base type of all binary components."""
abstract type BinaryComponent <: DelayComponent end

"""Abstract base type of all uncorrelated (white) noise components."""
abstract type WhiteNoiseComponent <: Component end

show(io::IO, ::MIME"text/plain", comp::Component) = show(io, comp)
