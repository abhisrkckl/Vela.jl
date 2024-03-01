using GeometricUnits

export Parameter,
    MultiParameter, read_param, ParamHandler, read_params, get_free_param_names

struct Parameter
    display_name::String
    default_quantity::GQ{Float64}
    frozen::Bool

    function Parameter(name, default_quantity, frozen)
        return new(uppercase(name), default_quantity, frozen)
    end
end

function read_param(param::Parameter, value::Float64)
    @assert !param.frozen "Refusing to read a frozen parameter $(param.display_name)."

    return quantity_like(param.default_quantity, value)
end

struct MultiParameter
    name::String
    parameters::Vector{Parameter}
end

struct ParamHandler
    multi_params::Vector{MultiParameter}
    _default_params_dict::Dict{String,Vector{GQ{Float64}}}

    function ParamHandler(multi_params)
        _default_params_dict = Dict{String,Vector{GQ{Float64}}}()

        for mpar in multi_params
            _default_params_dict[mpar.name] =
                [param.default_quantity for param in mpar.parameters]
        end

        return new(multi_params, _default_params_dict)
    end
end

function read_params(param_handler::ParamHandler, values::Vector{Float64})
    param_dict = deepcopy(param_handler._default_params_dict)
    ii = 1
    for mpar in param_handler.multi_params
        for (jj, param) in enumerate(mpar.parameters)
            if !param.frozen
                @inbounds param_dict[mpar.name][jj] = read_param(param, values[ii])
                ii += 1
            end
        end
    end

    return param_dict
end

function get_free_param_names(param_handler::ParamHandler)
    pnames = Vector{String}()
    ii = 1
    for mpar in param_handler.multi_params
        for param in mpar.parameters
            if !param.frozen
                @inbounds push!(pnames, param.display_name)
                ii += 1
            end
        end
    end

    return pnames
end
