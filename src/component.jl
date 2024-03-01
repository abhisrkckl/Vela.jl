using .Threads

export Component,
    UncorrelatedComponent,
    CorrelatedComponent,
    PhaseComponent,
    DelayComponent,
    UncorrelatedNoiseComponent,
    CorrelatedNoiseComponent,
    correct_toa

abstract type Component end

abstract type UncorrelatedComponent <: Component end
abstract type CorrelatedComponent <: Component end

abstract type PhaseComponent <: UncorrelatedComponent end

correct_toa(component::PhaseComponent, toa::TOA, params::NamedTuple) =
    correct_toa_phase(toa, phase(component, toa, params))

abstract type DelayComponent <: UncorrelatedComponent end

correct_toa(component::DelayComponent, toa::TOA, params::NamedTuple) =
    correct_toa_delay(toa, delay(component, toa, params))

abstract type DispersionComponent <: DelayComponent end

delay(component::DispersionComponent, toa::TOA, params) =
    dispersion_slope(component, toa, params) / toa.observing_frequency^2

# abstract type BasisDelayComponent <: DelayComponent end
# delay(component::BasisDelayComponent, toa::TOA, params) =
#     dot(basis(component, toa), amplitudes(component, toa, params))

# abstract type UncorrelatedNoiseComponent <: UncorrelatedComponent end

# abstract type CorrelatedNoiseComponent <: CorrelatedComponent end
