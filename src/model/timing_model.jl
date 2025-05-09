export TimingModel, correct_toa

"""
    TimingModel

The pulsar timing & noise model. Supports both narrowband and wideband timing.
Corresponds to `TimingModel` in `PINT`.

References:
    [Edwards+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10870.x),
    [Lentati+ 2014](http://doi.org/10.1093/mnras/stt2122),
    [Alam+ 2021](http://doi.org/10.3847/1538-4365/abc6a1)
"""
struct TimingModel{ComponentsTuple<:Tuple,KernelType<:Kernel,PriorsTuple<:Tuple}
    pulsar_name::String
    ephem::String
    clock::String
    units::String
    epoch::GQ{1,Float64}
    components::ComponentsTuple
    kernel::KernelType
    param_handler::ParamHandler
    tzr_toa::TOA
    priors::PriorsTuple

    function TimingModel(
        pulsar_name::String,
        ephem::String,
        clock::String,
        units::String,
        epoch::GQ{1,Float64},
        components::Tuple,
        kernel::Kernel,
        param_handler::ParamHandler,
        tzr_toa::TOA,
        priors::Tuple,
    )
        @assert is_tzr(tzr_toa) "Invalid TZR TOA."
        @assert all(map(c -> isa(c, Component), components)) "components must be a tuple containing Component objects."
        @assert all(map(p -> isa(p, Prior), priors)) "priors must be a tuple containing Prior objects."
        @assert length(priors) == length(get_free_param_names(param_handler)) "The number of prior distributions must match the number of free parameters."

        return new{typeof(components),typeof(kernel),typeof(priors)}(
            pulsar_name,
            ephem,
            clock,
            units,
            epoch,
            components,
            kernel,
            param_handler,
            tzr_toa,
            priors,
        )
    end
end

read_params(model::TimingModel, values) = read_params(model.param_handler, values)

get_free_param_names(model::TimingModel) = get_free_param_names(model.param_handler)
get_free_param_units(model::TimingModel) = get_free_param_units(model.param_handler)
get_free_param_labels(model::TimingModel) = get_free_param_labels(model.param_handler)
get_free_param_prefixes(model::TimingModel) = get_free_param_prefixes(model.param_handler)
get_num_timing_params(model::TimingModel) = get_num_timing_params(model.param_handler)
get_marginalized_param_names(model::TimingModel) =
    get_marginalized_param_names(model.kernel)

read_param_values_to_vector(model::TimingModel) =
    read_param_values_to_vector(model.param_handler)
read_param_values_to_vector(model::TimingModel, params::NamedTuple) =
    read_param_values_to_vector(model.param_handler, params)

get_scale_factors(model::TimingModel) = get_scale_factors(model.param_handler)

@unroll function correct_toa( # COV_EXCL_LINE
    components::Tuple,
    toa::TOABase,
    toacorr::TOACorrectionBase,
    params::NamedTuple,
)::TOACorrectionBase
    toacorr1 = toacorr
    @unroll for component in components
        toacorr1 = correct_toa(component, toa, toacorr1, params)
    end
    return toacorr1
end

"""Update a `CorrectedTOA` object using a timing model.
This involves applying the correction from each component in succession."""
correct_toa(
    model::TimingModel,
    toa::TOABase,
    toacorr::TOACorrectionBase,
    params::NamedTuple,
) = correct_toa(model.components, toa, toacorr, params)

"""Update a `TOA` object using a timing model."""
correct_toa(model::TimingModel, toa::TOA, params::NamedTuple) =
    correct_toa(model, toa, TOACorrection(), params)

show(io::IO, model::TimingModel) =
    print(io, "TimingModel[$(string(model.components)); $(string(model.kernel))]")
show(io::IO, ::MIME"text/plain", model::TimingModel) = show(io, model)
