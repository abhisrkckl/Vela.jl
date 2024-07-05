using GeometricUnits
using Unrolled
using .Threads

export TimingModel,
    correct_toa, form_residual, calc_chi2, calc_lnlike, calc_tzr_phase, get_lnlike_func

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

function form_residual(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase::GQ)::GQ
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
    return dphase / doppler_shifted_spin_frequency(ctoa)
end

function form_residuals(
    model::TimingModel,
    toas::Vector{TOA},
    params::NamedTuple,
)::Vector{GQ}
    tzrphase = calc_tzr_phase(model, params)
    return [form_residual(model, toa, params, tzrphase) for toa in toas]
end

function calc_tzr_phase(model::TimingModel, params::NamedTuple)
    ctzrtoa = correct_toa(model, model.tzr_toa, params)
    return ctzrtoa.phase
end

function calc_chi2(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    chisq = Atomic{Float64}(0.0)
    @threads for toa in toas
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(ctoa)

        atomic_add!(chisq, Float64((tres * tres / err2).x))
    end

    return chisq[]
end

function calc_chi2_serial(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    chisq = 0.0
    for toa in toas
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(ctoa)

        chisq += (tres * tres / err2).x
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
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(ctoa)
        norm = log(err2.x)

        atomic_add!(result, (tres * tres / err2).x + norm)
    end

    return -result[] / 2
end

function calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    result = 0.0
    for toa in toas
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(ctoa)
        norm = log(err2.x)

        result += (tres * tres / err2).x + norm
    end

    return -result / 2
end

calc_lnlike(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_lnlike(model, toas, read_params(model, params))

calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_lnlike_serial(model, toas, read_params(model, params))

show(io::IO, model::TimingModel) = print(io, "TimingModel[$(string(model.components))]")
show(io::IO, ::MIME"text/plain", model::TimingModel) = show(io, model)
