using GeometricUnits

export TimingModel

struct TimingModel
    components::Vector{Component}
    param_handler::ParamHandler
    tzr_toa::TOA

    function TimingModel(components::Components, param_handler::ParamHandler, tzr_toa::TOA)
        @assert tzr_toa.tzr "Invalid TZR TOA."

        return new(components, param_handler, tzr_toa)
    end
end
