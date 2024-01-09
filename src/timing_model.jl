using GeometricUnits

export TimingModel

struct TimingModel
    param_handler::ParamHandler
    tzr_toa::TOA

    function TimingModel(param_handler::ParamHandler, tzr_toa::TOA)
        @assert tzr_toa.tzr "Invalid TZR TOA."

        return new(param_handler, tzr_toa)
    end
end
