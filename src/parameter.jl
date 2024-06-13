using GeometricUnits

export Parameter,
    MultiParameter,
    read_param,
    ParamHandler,
    read_params,
    get_free_param_names,
    read_param_values_to_vector

struct Parameter
    name::Symbol
    default_quantity::GQ{Float64}
    frozen::Bool
    original_units::String
    unit_conversion_factor::Float64
end

struct MultiParameter
    name::Symbol
    parameters::Vector{Parameter}
end

function read_param(param::Parameter, value::Float64)
    @assert !param.frozen "Refusing to read a frozen parameter $(param.name)."

    return quantity_like(param.default_quantity, value)
end

function _get_single_params_tuple(single_params)
    pkeys = Tuple((sp.name for sp in single_params))
    pquants = Tuple((sp.default_quantity for sp in single_params))
    return (; zip(pkeys, pquants)...)
end

function _get_multi_params_tuple(multi_params)
    pkeys = Tuple((mp.name for mp in multi_params))
    ptups = Tuple((Tuple(_get_single_params_tuple(mp.parameters)) for mp in multi_params))
    return (; zip(pkeys, ptups)...)
end

_get_params_tuple(single_params, multi_params) =
    merge(_get_single_params_tuple(single_params), _get_multi_params_tuple(multi_params))

struct ParamHandler{ParamsType<:NamedTuple}
    single_params::Vector{Parameter}
    multi_params::Vector{MultiParameter}
    _default_param_quantities::ParamsType
    _default_quantities::Vector{GQ{Float64}}
    _free_indices::Vector{Int}
end

function ParamHandler(single_params, multi_params)
    default_params = _get_params_tuple(single_params, multi_params)
    default_quantities = collect(Iterators.flatten(default_params))

    all_params = [
        single_params
        collect(Iterators.flatten([mpar.parameters for mpar in multi_params]))
    ]
    free_indices = collect(1:length(all_params))[[!par.frozen for par in all_params]]

    ParamHandler(
        single_params,
        multi_params,
        default_params,
        default_quantities,
        free_indices,
    )
end

function read_params(
    ph::ParamHandler{ParamsType},
    free_values::Vector{Float64},
)::ParamsType where {ParamsType<:NamedTuple}
    quantities = copy(ph._default_quantities)
    for (idx, value) in zip(ph._free_indices, free_values)
        d = quantities[idx].d
        quantities[idx] = GQ(value, d)
    end

    return reinterpret(ParamsType, quantities)[1]
end

# function read_params(param_handler::ParamHandler, values::Vector{Float64})
#     param_dict = deepcopy(param_handler._default_params_dict)
#     ii = 1
#     for mpar in param_handler.multi_params
#         for (jj, param) in enumerate(mpar.parameters)
#             if !param.frozen
#                 @inbounds param_dict[mpar.name][jj] = read_param(param, values[ii])
#                 ii += 1
#             end
#         end
#     end

#     return param_dict
# end

function get_free_param_names(param_handler::ParamHandler)
    pnames = Vector{String}()
    ii = 1
    for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                @inbounds push!(pnames, string(param.name))
                ii += 1
            end
        end
    end

    return pnames
end

function read_param_values_to_vector(param_handler::ParamHandler, params::NamedTuple)
    param_vec = Float64[]
    @inbounds for mpar in param_handler.multi_params
        if length(mpar.parameters) == 1
            if !mpar.parameters[1].frozen
                push!(param_vec, Float64(params[mpar.name].x))
            end
        else
            for (jj, param) in enumerate(mpar.parameters)
                if !param.frozen
                    push!(param_vec, Float64(params[mpar.name][jj].x))
                end
            end
        end

    end

    return param_vec
end
