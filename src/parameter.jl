using GeometricUnits

export Parameter, read_param, read_params

struct Parameter
    name::String
    default_quantity::GQ{Float64}
    upper_limit::GQ{Float64}
    lower_limit::GQ{Float64}
    frozen::Bool

    function Parameter(name, default_quantity, lower_limit, upper_limit, frozen)
        @assert default_quantity.d == upper_limit.d == lower_limit.d "upper_limit and lower_limit must have the same dimension as default_quantity."
        @assert lower_limit <= default_quantity <= upper_limit "default_quantity outside limits."

        return new(uppercase(name), default_quantity, upper_limit, lower_limit, frozen)
    end
end

function read_param(param::Parameter, value::Float64)
    @assert !param.frozen "Refusing to read a frozen parameter $(param.name)."

    q = quantity_like(param.default_quantity, value)
    @assert param.lower_limit <= q <= param.upper_limit "Given value is outside the limits for $(param.name)."

    return q
end

struct ParamHandler
    params::Vector{Parameter}
    _free_params::Vector{String}
    _default_params_dict::Dict{String,GQ{Float64}}

    function ParamHandler(params)
        param_names = [param.name for param in params]
        @assert length(param_names) == length(Set(param_names)) "Repeated parameter names found in default_params."

        _free_params = [param for param in params if !param.frozen]
        _default_params_dict =
            Dict(param.name => param.default_quantity for param in params)

        return new(default_params, _free_params, _default_params_dict)
    end
end

function read_params(param_handler::ParamHandler, values::Vector{Float64})
    free_params_dict = Dict(
        param.name => read_param(param, value) for
        (param, value) in zip(param_handler._free_params, values)
    )
    return merge(param_handler._default_params_dict, free_params_dict)
end
