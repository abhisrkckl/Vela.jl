using GeometricUnits
using .Threads

export TimingModel, correct_toa, form_residual, calc_chi2, calc_lnlike, calc_tzr_phase

struct TimingModel{ComponentsTuple<:Tuple}
    components::ComponentsTuple
    param_handler::ParamHandler
    tzr_toa::TOA

    function TimingModel(components::Tuple, param_handler::ParamHandler, tzr_toa::TOA)
        @assert tzr_toa.tzr "Invalid TZR TOA."
        @assert all(map(c -> isa(c, Component), components)) "components must be a tuple containing Component objects."

        return new{typeof(components)}(components, param_handler, tzr_toa)
    end
end

read_params_from_dict(model::TimingModel, params::Dict)::NamedTuple =
    merge((read_params_from_dict(component, params) for component in model.components)...)

read_params(model::TimingModel, values::Vector{Float64}) =
    return read_params_from_dict(model, read_params(model.param_handler, values))

function correct_toa(model::TimingModel, toa::TOA, params::NamedTuple)
    ctoa::TOA = toa
    for component in model.components
        ctoa = correct_toa(component, ctoa, params)
    end
    return ctoa
end

function form_residual(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase::GQ)::GQ
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(ctoa.phase - tzrphase)
    return dphase / ctoa.spin_frequency
end

function calc_tzr_phase(model::TimingModel, params::NamedTuple)
    ctzrtoa = correct_toa(model, model.tzr_toa, params)
    return ctzrtoa.phase
end

function calc_chi2(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    chisq = Atomic{Float64}(0.0)
    @threads for toa in toas
        res = form_residual(model, toa, params, tzrphase)
        err = toa.error

        atomic_add!(chisq, Float64(((res / err)^2).x))
    end

    return chisq[]
end

function calc_chi2_serial(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    chisq = 0.0
    for toa in toas
        res = form_residual(model, toa, params, tzrphase)
        err = toa.error

        chisq += ((res / err)^2).x
    end

    return chisq
end

calc_chi2(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_chi2(model, toas, read_params(model, params))

calc_chi2_serial(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_chi2_serial(model, toas, read_params(model, params))

function calc_lnlike(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    result = Atomic{Float64}(0.0)
    @threads for toa in toas
        res = form_residual(model, toa, params, tzrphase)
        err = toa.error
        norm = log(toa.error.x)

        atomic_add!(result, ((res / err)^2).x + norm)
    end

    return -result[] / 2
end

function calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    result = 0.0
    for toa in toas
        res = form_residual(model, toa, params, tzrphase)
        err = toa.error
        norm = log(toa.error.x)

        result += ((res / err)^2).x + norm
    end

    return -result / 2
end

calc_lnlike(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_lnlike(model, toas, read_params(model, params))

calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_lnlike_serial(model, toas, read_params(model, params))

show(io::IO, model::TimingModel) = print(io, "TimingModel[$(string(model.components))]")
show(io::IO, ::MIME"text/plain", model::TimingModel) = show(io, model)
