using GeometricUnits

export TimingModel, correct_toas, form_residuals, calc_chi2, calc_lnlike

struct TimingModel
    components::Vector{Component}
    param_handler::ParamHandler
    tzr_toa::TOA

    function TimingModel(
        components::Vector{Component},
        param_handler::ParamHandler,
        tzr_toa::TOA,
    )
        @assert tzr_toa.tzr "Invalid TZR TOA."

        return new(components, param_handler, tzr_toa)
    end
end

function correct_toas(model::TimingModel, toas::Vector{TOA}, params::Dict)
    corrected_toas::Vector{TOA} = copy(toas)

    for component in model.components
        correct_toas!(component, corrected_toas, params)
    end

    return corrected_toas
end

function form_residuals(model::TimingModel, toas::Vector{TOA}, params)
    corrected_tzr_toa = correct_toas(model, [model.tzr_toa], params)[1]
    corrected_toas = correct_toas(model, toas, params)
    return [
        (corrected_toa.phase - corrected_tzr_toa.phase) / corrected_toa.spin_frequency for
        corrected_toa in corrected_toas
    ]
end

function calc_chi2(model::TimingModel, toas::Vector{TOA}, params)
    resids = form_residuals(model, toas, params)
    return sum((resid / toa.error)^2 for (resid, toa) in zip(resids, toas))
end

function calc_lnlike(model::TimingModel, toas::Vector{TOA}, params)
    resids = form_residuals(model, toas, params)
    return -0.5 * sum(
        (resid / toa.error)^2 + log(toa.error.x) for (resid, toa) in zip(resids, toas)
    )
end
