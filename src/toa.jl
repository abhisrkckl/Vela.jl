using GeometricUnits
using Quadmath

export EphemerisVectors, TOA, is_barycentered, correct_toa_delay, correct_toa_phase, make_tzr_toa

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
    tzr::Bool
end

TOA(value, error, frequency, phase, ephem_vectors) =
    TOA(value, error, frequency, phase, ephem_vectors, 0, false)
TOA(value, error, frequency, phase) = TOA(value, error, frequency, phase, nothing, 0, false)

is_barycentered(toa::TOA) = isnothing(toa.ephem_vectors)

make_tzr_toa(tzrtdb, tzrfreq, ephem_vectors) = TOA(tzrtdb, time(0.0), tzrfreq, dimensionless(Float128(0.0)), ephem_vectors, 0, true)

correct_toa_delay(toa::TOA, delay::GQ) = TOA(
    toa.value - delay,
    toa.error,
    toa.frequency,
    toa.phase,
    toa.ephem_vectors,
    toa.level + 1,
    toa.tzr,
)
correct_toa_phase(toa::TOA, phase::GQ) = TOA(
    toa.value,
    toa.error,
    toa.frequency,
    toa.phase + phase,
    toa.ephem_vectors,
    toa.level + 1,
    toa.tzr,
)
