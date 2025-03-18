export Component,
    PhaseComponent,
    DelayComponent,
    WhiteNoiseComponent,
    BinaryComponent,
    DispersionComponent,
    correct_toa,
    delay,
    is_gp_noise

"""
    Component

Abstract base type of all timing & noise model components which affect one TOA at a time."""
abstract type Component end

"""
    correct_toa(::Component, ::TOABase, ::TOACorrectionBase, ::NamedTuple)

Correct the TOA using a delay, phase, observing frequency shift, uncertainty scaling, doppler factor, etc."""
function correct_toa end

"""
    is_gp_noise(::Component)::Bool

Whether a component represents a correlated Gaussian noise process.
"""
is_gp_noise(::Component) = false

"""
    PhaseComponent

Abstract base type of all timing model components which contribute a phase correction to a TOA."""
abstract type PhaseComponent <: Component end

correct_toa(
    component::PhaseComponent,
    toa::TOABase,
    toacorr::TOACorrection,
    params::NamedTuple,
) = correct_toa_phase(toacorr; phase = phase(component, toa, toacorr, params))

"""
    DelayComponent

Abstract base type of all timing model components which contribute a time delay correction to a TOA."""
abstract type DelayComponent <: Component end

correct_toa(
    component::DelayComponent,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
) = correct_toa_delay(toacorr; delay = delay(component, toa, toacorr, params))

"""
    DispersionComponent

Abstrct base type of all timing model components which provide a dispersion measure correction."""
abstract type DispersionComponent <: DelayComponent end

"""
    dispersion_slope(::DispersionComponent, ::TOA, ::TOACorrection, ::NamedTuple)

Compute the dispersion slope corresponding to a TOA. 
"""
function dispersion_slope end

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

"""
    BinaryComponent

Abstract base type of all binary components."""
abstract type BinaryComponent <: DelayComponent end

"""
    WhiteNoiseComponent

Abstract base type of all uncorrelated (white) noise components."""
abstract type WhiteNoiseComponent <: Component end

show(io::IO, ::MIME"text/plain", comp::Component) = show(io, comp)
