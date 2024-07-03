export correct_toa, form_residual

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
