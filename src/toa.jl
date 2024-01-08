using GeometricUnits
using Quadmath
using LinearAlgebra

export EphemerisVectors,
    TOA, is_barycentered, correct_toa_delay, correct_toa_phase, make_tzr_toa

struct EphemerisVectors
    ssb_obs_pos::Vector{GQ{Float64}}
    ssb_obs_vel::Vector{GQ{Float64}}
    obs_sun_pos::Vector{GQ{Float64}}

    function EphemerisVectors(ssb_obs_pos, ssb_obs_vel, obs_sun_pos)
        @assert all([xi.d == 1 for xi in ssb_obs_pos]) "Dimension mismatch in ssb_obs_pos."
        @assert all([vi.d == 0 for vi in ssb_obs_vel]) "Dimension mismatch in ssb_obs_vel."
        @assert all([vi.d == 1 for vi in obs_sun_pos]) "Dimension mismatch in obs_sun_pos."
        @assert dot(ssb_obs_vel, ssb_obs_vel) < 1 "Magnitude of ssb_obs_vel should be less than 1."

        return new(ssb_obs_pos, ssb_obs_vel, obs_sun_pos)
    end
end

struct TOA
    value::GQ{Float128}
    error::GQ{Float64}
    frequency::GQ{Float64}
    phase::GQ{Float128}
    ephem_vectors::Union{EphemerisVectors,Nothing}
    level::UInt
    tzr::Bool

    function TOA(value, error, frequency, phase, ephem_vectors, level, tzr)
        @assert value.d == 1 "Dimension mismatch in value (given $(value.d), expected 1)."
        @assert error.d == 1 "Dimension mismatch in error (given $(error.d), expected 1)."
        @assert frequency.d == -1 "Dimension mismatch in frequency (given $(frequency.d), expected -1)."
        @assert phase.d == 0 "Dimension mismatch in phase (given $(phase.d), expected 0)."

        return new(value, error, frequency, phase, ephem_vectors, level, tzr)
    end
end

TOA(value, error, frequency, phase, ephem_vectors) =
    TOA(value, error, frequency, phase, ephem_vectors, 0, false)
TOA(value, error, frequency, phase) = TOA(value, error, frequency, phase, nothing, 0, false)

is_barycentered(toa::TOA) = isnothing(toa.ephem_vectors)

make_tzr_toa(tzrtdb, tzrfreq, ephem_vectors) =
    TOA(tzrtdb, time(0.0), tzrfreq, dimensionless(Float128(0.0)), ephem_vectors, 0, true)

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
