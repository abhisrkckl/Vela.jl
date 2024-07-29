export TimingModel, correct_toa, form_residual, calc_chi2, calc_lnlike, calc_tzr_phase

"""The pulsar timing & noise model.

Corresponds to `TimingModel` in `PINT`."""
struct TimingModel{ComponentsTuple<:Tuple,PriorsTuple<:Tuple}
    pulsar_name::String
    ephem::String
    clock::String
    units::String
    components::ComponentsTuple
    param_handler::ParamHandler
    tzr_toa::TOA
    priors::PriorsTuple

    function TimingModel(
        pulsar_name::String,
        ephem::String,
        clock::String,
        units::String,
        components::Tuple,
        param_handler::ParamHandler,
        tzr_toa::TOA,
        priors::Tuple,
    )
        @assert tzr_toa.tzr "Invalid TZR TOA."
        @assert all(map(c -> isa(c, Component), components)) "components must be a tuple containing Component objects."
        @assert all(map(p -> isa(p, Prior), priors)) "priors must be a tuple containing Prior objects."
        @assert length(priors) == length(get_free_param_names(param_handler)) "The number of prior distributions must match the number of free parameters."

        return new{typeof(components),typeof(priors)}(
            pulsar_name,
            ephem,
            clock,
            units,
            components,
            param_handler,
            tzr_toa,
            priors,
        )
    end
end

"""Create a parameter tuple from a collection of free parameter values."""
read_params(model::TimingModel, values) = read_params(model.param_handler, values)

@unroll function correct_toa(
    components::Tuple,
    ctoa::CorrectedTOA,
    params::NamedTuple,
)::CorrectedTOA
    ctoa1 = ctoa
    @unroll for component in components
        ctoa1 = correct_toa(component, ctoa1, params)
    end
    return ctoa1
end

"""Update a `CorrectedTOA` object using a timing model.
This involves applying the correction from each component in succession."""
correct_toa(model::TimingModel, ctoa::CorrectedTOA, params::NamedTuple) =
    correct_toa(model.components, ctoa, params)

"""Update a `TOA` object using a timing model."""
correct_toa(model::TimingModel, toa::TOA, params::NamedTuple) =
    correct_toa(model, CorrectedTOA(toa), params)

show(io::IO, model::TimingModel) = print(io, "TimingModel[$(string(model.components))]")
show(io::IO, ::MIME"text/plain", model::TimingModel) = show(io, model)
