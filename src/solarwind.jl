export SolarWind, SolarWindDispersion

abstract type SolarWind <: DispersionComponent end

struct SolarWindDispersion <: SolarWind
    model::Int

    function SolarWindDispersion(model::Int)
        @assert (model in [0, 1]) "Invalid solar wind model!"
        return new(model)
    end
end

dispersion_slope(swd::SolarWindDispersion, toa::TOA, params) = GQ(0.0, -1)
