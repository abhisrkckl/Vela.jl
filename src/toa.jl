using GeometricUnits
using Quadmath

export TOA, is_barycentered, correct_toa_delay, correct_toa_phase

struct EphemerisVectors
    ssb_obs_pos::Vector{GQ{Float64}}
    ssb_obs_vel::Vector{GQ{Float64}}
    obs_sun_pos::Vector{GQ{Float64}}
end

struct TOA
    value::GQ{Float128}
    error::GQ{Float64}
    frequency::GQ{Float64}
    phase::GQ{Float128}
    ephem_vectors::Union{EphemerisVectors,Nothing}
    level::UInt
end

TOA(value, error, frequency, phase, ephem_vectors) = TOA(value, error, frequency, phase, ephem_vectors, 0)
TOA(value, error, frequency, phase) = TOA(value, error, frequency, phase, nothing, 0)

is_barycentered(toa::TOA) = isnothing(toa.ephem_vectors)

correct_toa_delay(toa::TOA, delay::GQ) = TOA(toa.value - delay, toa.error, toa.frequency, toa.phase, toa.ephem_vectors, toa.level+1)
correct_toa_phase(toa::TOA, phase::GQ) = TOA(toa.value, toa.error, toa.frequency, toa.phase + phase, toa.ephem_vectors, toa.level+1)