export SolarSystem

struct SolarSystem <: DelayComponent
    ecliptic_coordinates::Bool
    planet_shapiro::Bool
end

delay(::SolarSystem, toa::TOA, params) = 0
