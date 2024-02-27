using .Threads

export Component,
    UncorrelatedComponent,
    CorrelatedComponent,
    PhaseComponent,
    DelayComponent,
    UncorrelatedNoiseComponent,
    CorrelatedNoiseComponent,
    correct_toas!

abstract type Component end

abstract type UncorrelatedComponent <: Component end
abstract type CorrelatedComponent <: Component end

abstract type PhaseComponent <: UncorrelatedComponent end

function correct_toas!(component::PhaseComponent, toas::Vector{TOA}, params::Dict)
    params_tuple = read_params_from_dict(component, params)

    @threads for toa in toas
        correct_toa_phase!(toa, phase(component, toa, params_tuple))
    end
end

abstract type DelayComponent <: UncorrelatedComponent end

function correct_toas!(component::DelayComponent, toas::Vector{TOA}, params::Dict)
    params_tuple = read_params_from_dict(component, params)

    @threads for toa in toas
        correct_toa_delay!(toa, delay(component, toa, params_tuple))
    end
end

abstract type DispersionComponent <: DelayComponent end

delay(component::DispersionComponent, toa::TOA, params) =
    dispersion_slope(component, toa, params) / toa.observing_frequency^2

# abstract type BasisDelayComponent <: DelayComponent end
# delay(component::BasisDelayComponent, toa::TOA, params) =
#     dot(basis(component, toa), amplitudes(component, toa, params))

# abstract type UncorrelatedNoiseComponent <: UncorrelatedComponent end

# abstract type CorrelatedNoiseComponent <: CorrelatedComponent end
