export Component,
    UncorrelatedComponent,
    CorrelatedComponent,
    PhaseComponent,
    DelayComponent,
    UncorrelatedNoiseComponent,
    CorrelatedNoiseComponent,
    correct_toa

"""Abstract base type of all timing & noise model components."""
abstract type Component end

"""Abstract base type of all timing & noise model components which affect one TOA at a time."""
abstract type UncorrelatedComponent <: Component end

"""Abstract base type of all timing & noise model components which affect multiple TOAs at a time."""
abstract type CorrelatedComponent <: Component end

"""Abstract base type of all timing model components which contribute a phase correction to a TOA."""
abstract type PhaseComponent <: UncorrelatedComponent end

"""Correct the phase of a `TOA`."""
correct_toa(component::PhaseComponent, ctoa::CorrectedTOA, params::NamedTuple) =
    correct_toa(ctoa; phase = phase(component, ctoa, params))

"""Abstract base type of all timing model components which contribute a time delay correction to a TOA."""
abstract type DelayComponent <: UncorrelatedComponent end

"""Correct the value of a `TOA` using a delay."""
correct_toa(component::DelayComponent, ctoa::CorrectedTOA, params::NamedTuple) =
    correct_toa(ctoa; delay = delay(component, ctoa, params))

"""Abstract base type of all timing model components which contribute a time delay and a 
Doppler factor correction to a TOA."""
abstract type KinematicDelayComponent <: DelayComponent end

"""Abstrct base type of all timing model components which provide a dispersion measure correction."""
abstract type DispersionComponent <: DelayComponent end

"""Compute a dispersion delay."""
delay(component::DispersionComponent, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    dispersion_slope(component, ctoa, params) /
    doppler_corrected_observing_frequency(ctoa)^Val(2)

"""Abstract base type of all binary components."""
abstract type BinaryComponent <: KinematicDelayComponent end

"""Abstract base type of all uncorrelated (white) noise components."""
abstract type WhiteNoiseComponent <: UncorrelatedComponent end

show(io::IO, ::MIME"text/plain", comp::Component) = show(io, comp)
