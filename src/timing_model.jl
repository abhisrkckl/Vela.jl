export TimingModel

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

function calc_lnprior(model::TimingModel, params::Vector{Float64})
    result = 0
    for component in model.components
        priors = prior_distributions(component)
        result += logpdf(priors, params)
    end
    return result
end

show(io::IO, model::TimingModel) = print(io, "TimingModel[$(string(model.components))]")
show(io::IO, ::MIME"text/plain", model::TimingModel) = show(io, model)
