using GeometricUnits

export TimingModel

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

function correct_toa(model::TimingModel, toa::TOA, params)::TOA
    corrected_toa = toa

    for component in model.components
        corrected_toa = correct_toa(component, corrected_toa, params)
    end

    return corrected_toa
end

function form_residuals(model::TimingModel, toas::Vector{TOA}, params)
    corrected_tzr_toa = correct_toa(model, model.tzr_toa, params)
    corrected_toas = [correct_toa(model, toa, params) for toa in toas]
    return [
        (corrected_toa.phase - corrected_tzr_toa.phase) / correct_toa.spin_frequency for
        corrected_toa in corrected_toas
    ]
end

function calc_chi2(model::TimingModel, toas::Vector{TOA}, params)
    resids = form_residuals(model, toas, params)
    return sum((resid / toa.error)^2 for (resid, toa) in zip(resids, toas))
end
