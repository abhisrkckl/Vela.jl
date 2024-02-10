export SolarSystem, delay

struct SolarSystem <: DelayComponent
    ecliptic_coordinates::Bool
    planet_shapiro::Bool
end

delay(::SolarSystem, toa::TOA, params) = time(0.0)