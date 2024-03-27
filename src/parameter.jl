using GeometricUnits
using Accessors

export Parameter,
    MultiParameter,
    read_param,
    ParamHandler,
    read_params,
    get_free_param_names,
    read_param_values_to_vector

struct Parameter
    display_name::Symbol
    default_quantity::GQ{Float64}
    frozen::Bool

    function Parameter(name, default_quantity, frozen)
        return new(name, default_quantity, frozen)
    end
end

function read_param(param::Parameter, value::Float64)
    @assert !param.frozen "Refusing to read a frozen parameter $(param.display_name)."

    return quantity_like(param.default_quantity, value)
end

struct MultiParameter
    name::Symbol
    parameters::Vector{Parameter}
end

struct ParamHandler
    multi_params::Vector{MultiParameter}
    _params_tuple::NamedTuple

    function ParamHandler(multi_params)
        # _default_quantities = collect(Iterators.flatten(map(mp -> map(p -> p.default_quantity, mp.parameters), multi_params)))
        
        keys = Tuple([mpar.name for mpar in multi_params])
        types = [
            (
                length(mpar.parameters) == 1 ? 
                typeof(mpar.parameters[1].default_quantity) : 
                NTuple{length(mpar.parameters), typeof(mpar.parameters[1].default_quantity)}
            )
            for mpar in multi_params
        ]
        args = [
            (
                length(mpar.parameters) == 1 ? 
                mpar.parameters[1].default_quantity : 
                Tuple([param.default_quantity for param in mpar.parameters])
            )
            for mpar in multi_params
        ]

        _params_tuple = NamedTuple{keys, Tuple{types...,}}(args)

        return new(multi_params, _params_tuple)
    end
end

function read_params(param_handler::ParamHandler, values::Vector{Float64})
    result = param_handler._params_tuple
    jj = 1
    for (ii, mpar) in enumerate(param_handler.multi_params)
        if length(mpar.parameters) == 1 
            if !mpar.parameters[1].frozen
                @reset result[ii] = quantity_like(result[ii], values[jj])
                jj += 1
            end
        else 
            for (ij, param) in enumerate(mpar.parameters)
                if !param.frozen
                    @reset result[ii][ij] = quantity_like(result[ii][ij], values[jj])
                    jj += 1
                end
            end
        end
    end

    return result
end

function get_free_param_names(param_handler::ParamHandler)
    pnames = Vector{String}()
    ii = 1
    for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                @inbounds push!(pnames, string(param.display_name))
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