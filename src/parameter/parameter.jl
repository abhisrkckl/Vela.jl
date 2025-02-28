export Parameter,
    MultiParameter,
    read_param,
    ParamHandler,
    read_params,
    get_free_param_names,
    get_free_param_units,
    get_free_param_labels,
    get_free_param_prefixes,
    read_param_values_to_vector,
    get_scale_factors,
    is_noise

"""
    Parameter

A single model parameter.

Corresponds to `floatParameter`, `AngleParameter`, or `MJDParameter` in `PINT`.
"""
struct Parameter
    name::Symbol
    default_value::Float64
    dimension::Int
    frozen::Bool
    original_units::String
    unit_conversion_factor::Float64
    noise::Bool
end

Parameter(
    name::Symbol,
    default_quantity::GQ,
    frozen::Bool,
    original_units::String,
    unit_conversion_factor::Float64,
    noise::Bool,
) = Parameter(
    name,
    value(default_quantity),
    udim(default_quantity),
    frozen,
    original_units,
    unit_conversion_factor,
    noise,
)

"""
    MultiParameter

A set of model parameters that are characterized by a common name and a varying index.

Corresponds to `maskParameter` or `prefixParameter` in `PINT`.
"""
struct MultiParameter
    name::Symbol
    parameters::Vector{Parameter}

    function MultiParameter(name::Symbol, parameters::Vector{Parameter})
        @assert allequal([par.noise for par in parameters])
        return new(name, parameters)
    end
end

is_noise(mpar::MultiParameter) = mpar.parameters[1].noise

function _get_single_params_tuple(single_params)
    pkeys = Tuple((sp.name for sp in single_params))
    pquants = Tuple((GQ{sp.dimension}(sp.default_value) for sp in single_params))
    return (; zip(pkeys, pquants)...)
end

function _get_multi_params_tuple(multi_params)
    pkeys = Tuple((mp.name for mp in multi_params))
    ptups = Tuple((Tuple(_get_single_params_tuple(mp.parameters)) for mp in multi_params))
    return (; zip(pkeys, ptups)...)
end

_get_params_tuple(single_params, multi_params) =
    merge(_get_single_params_tuple(single_params), _get_multi_params_tuple(multi_params))

"""
    ParamHandler{PT<:NamedTuple}

Handles the creation of a parameter tuple from a collection of free parameter values.
"""
struct ParamHandler{ParamsType<:NamedTuple}
    single_params::Vector{Parameter}
    multi_params::Vector{MultiParameter}
    _default_params_tuple::ParamsType
    _default_values::Vector{Float64}
    _free_indices::Vector{Int}
    _nfree::Int
end

function ParamHandler(single_params, multi_params)
    default_params = _get_params_tuple(single_params, multi_params)
    default_values = collect(map(value, Iterators.flatten(default_params)))

    all_params = [
        single_params
        collect(Iterators.flatten([mpar.parameters for mpar in multi_params]))
    ]
    free_indices = collect(1:length(all_params))[[!par.frozen for par in all_params]]

    ParamHandler(
        single_params,
        multi_params,
        default_params,
        default_values,
        free_indices,
        length(free_indices),
    )
end

"""
    read_params(ph::ParamHandler{PT}, free_values)::PT where {PT<:NamedTuple}

Create a parameter tuple from a collection of free parameter values.

Reverse of `read_param_values_to_vector()`.
"""
function read_params(
    ph::ParamHandler{ParamsType},
    free_values,
)::ParamsType where {ParamsType<:NamedTuple}
    @assert length(free_values) == ph._nfree
    values = copy(ph._default_values)
    @inbounds for (idx, value) in zip(ph._free_indices, free_values)
        values[idx] = value
    end

    return reinterpret(ParamsType, values)[1]
end

"""Generate an ordered collection of free parameter names."""
function get_free_param_names(param_handler::ParamHandler)::Vector{String}
    pnames = Vector{String}()

    @inbounds for spar in param_handler.single_params
        if !spar.frozen
            push!(pnames, string(spar.name))
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                push!(pnames, string(param.name))
            end
        end
    end

    return pnames
end


"""Generate an ordered collection of free parameter prefixes."""
function get_free_param_prefixes(param_handler::ParamHandler)::Vector{String}
    pnames = Vector{String}()

    @inbounds for spar in param_handler.single_params
        if !spar.frozen
            push!(pnames, string(spar.name))
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                push!(pnames, string(mpar.name))
            end
        end
    end

    return pnames
end

"""Generate an ordered collection of free parameter labels."""
function get_free_param_labels(
    param_handler::ParamHandler;
    delim::String = "\n",
)::Vector{String}
    pnames = Vector{String}()

    @inbounds for spar in param_handler.single_params
        if !spar.frozen
            push!(
                pnames,
                isempty(spar.original_units) ? string(spar.name) :
                "$(string(spar.name))$delim($(spar.original_units))",
            )
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                push!(
                    pnames,
                    isempty(param.original_units) ? string(param.name) :
                    "$(string(param.name))$delim($(param.original_units))",
                )
            end
        end
    end

    return pnames
end

"""Generate an ordered collection of free parameter units."""
function get_free_param_units(param_handler::ParamHandler)::Vector{String}
    pnames = Vector{String}()

    @inbounds for spar in param_handler.single_params
        if !spar.frozen
            push!(pnames, spar.original_units)
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                push!(pnames, param.original_units)
            end
        end
    end

    return pnames
end

"""Generate a collection of free parameter values from a parameter tuple.

Reverse of `read_params()`
"""
function read_param_values_to_vector(
    param_handler::ParamHandler,
    params::NamedTuple,
)::Vector{Float64}
    param_vec = Float64[]

    @inbounds for spar in param_handler.single_params
        if !spar.frozen
            push!(param_vec, Float64(params[spar.name].x))
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for (jj, param) in enumerate(mpar.parameters)
            if !param.frozen
                push!(param_vec, Float64(params[mpar.name][jj].x))
            end
        end
    end

    return param_vec
end
read_param_values_to_vector(param_handler::ParamHandler) =
    read_param_values_to_vector(param_handler, param_handler._default_params_tuple)

"""Get the scale factors that convert the free parameters from `Vela.jl`'s 
internal representation to the units used in `PINT`."""
function get_scale_factors(param_handler::ParamHandler)
    scale_factors = Float64[]

    @inbounds for spar in param_handler.single_params
        if !spar.frozen
            push!(scale_factors, spar.unit_conversion_factor)
        end
    end

    @inbounds for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                push!(scale_factors, param.unit_conversion_factor)
            end
        end
    end

    return scale_factors
end
