export TimingModel, correct_toa, form_residual, calc_chi2, calc_lnlike, calc_tzr_phase

struct TimingModel{ComponentsTuple<:Tuple}
    pulsar_name::String
    ephem::String
    clock::String
    units::String
    components::ComponentsTuple
    param_handler::ParamHandler
    tzr_toa::TOA

    function TimingModel(
        pulsar_name::String,
        ephem::String,
        clock::String,
        units::String,
        components::Tuple,
        param_handler::ParamHandler,
        tzr_toa::TOA,
    )
        @assert tzr_toa.tzr "Invalid TZR TOA."
        @assert all(map(c -> isa(c, Component), components)) "components must be a tuple containing Component objects."

        return new{typeof(components)}(
            pulsar_name,
            ephem,
            clock,
            units,
            components,
            param_handler,
            tzr_toa,
        )
    end
end

read_params(model::TimingModel, values::Vector{Float64}) =
    return read_params(model.param_handler, values)

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

correct_toa(model::TimingModel, ctoa::CorrectedTOA, params::NamedTuple) =
    correct_toa(model.components, ctoa, params)

correct_toa(model::TimingModel, toa::TOA, params::NamedTuple) =
    correct_toa(model, CorrectedTOA(toa), params)

show(io::IO, model::TimingModel) = print(io, "TimingModel[$(string(model.components))]")
show(io::IO, ::MIME"text/plain", model::TimingModel) = show(io, model)
