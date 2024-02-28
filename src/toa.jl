using GeometricUnits
using Quadmath
using LinearAlgebra
import Base.copy

export EphemerisVectors,
    TOA, is_barycentered, correct_toa_delay!, correct_toa_phase!, make_tzr_toa

struct EphemerisVectors
    ssb_obs_pos::Vector{GQ{Float64}}
    ssb_obs_vel::Vector{GQ{Float64}}
    obs_sun_pos::Vector{GQ{Float64}}
    obs_jupiter_pos::Vector{GQ{Float64}}
    obs_saturn_pos::Vector{GQ{Float64}}
    obs_venus_pos::Vector{GQ{Float64}}
    obs_uranus_pos::Vector{GQ{Float64}}
    obs_neptune_pos::Vector{GQ{Float64}}
    obs_earth_pos::Vector{GQ{Float64}}

    function EphemerisVectors(
        ssb_obs_pos,
        ssb_obs_vel,
        obs_sun_pos,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
        obs_earth_pos,
    )
        @assert all([xi.d == 1 for xi in ssb_obs_pos]) "Dimension mismatch in ssb_obs_pos."
        @assert all([vi.d == 0 for vi in ssb_obs_vel]) "Dimension mismatch in ssb_obs_vel."
        @assert all([xi.d == 1 for xi in obs_sun_pos]) "Dimension mismatch in obs_sun_pos."
        @assert all([xi.d == 1 for xi in obs_jupiter_pos]) "Dimension mismatch in obs_jupiter_pos."
        @assert all([xi.d == 1 for xi in obs_saturn_pos]) "Dimension mismatch in obs_saturn_pos."
        @assert all([xi.d == 1 for xi in obs_venus_pos]) "Dimension mismatch in obs_venus_pos."
        @assert all([xi.d == 1 for xi in obs_uranus_pos]) "Dimension mismatch in obs_uranus_pos."
        @assert all([xi.d == 1 for xi in obs_neptune_pos]) "Dimension mismatch in obs_neptune_pos."
        @assert all([xi.d == 1 for xi in obs_earth_pos]) "Dimension mismatch in obs_earth_pos."
        @assert dot(ssb_obs_vel, ssb_obs_vel) < 1 "Magnitude of ssb_obs_vel should be less than 1."

        return new(
            ssb_obs_pos,
            ssb_obs_vel,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )
    end
end

mutable struct TOA
    value::GQ{Float128}
    error::GQ{Float64}
    observing_frequency::GQ{Float64}
    phase::GQ{Float128}
    spin_frequency::GQ{Float64}
    ephem_vectors::EphemerisVectors
    level::UInt
    tzr::Bool

    function TOA(
        value,
        error,
        observing_frequency,
        phase,
        spin_frequency,
        ephem_vectors,
        level,
        tzr,
    )
        @assert value.d == 1 "Dimension mismatch in value (given $(value.d), expected 1)."
        @assert error.d == 1 "Dimension mismatch in error (given $(error.d), expected 1)."
        @assert observing_frequency.d == -1 "Dimension mismatch in observing_frequency (given $(observing_frequency.d), expected -1)."
        @assert isnothing(spin_frequency) || spin_frequency.d == -1 "Dimension mismatch in spin_frequency (given $(spin_frequency.d), expected -1)."
        @assert phase.d == 0 "Dimension mismatch in phase (given $(phase.d), expected 0)."

        return new(
            value,
            error,
            observing_frequency,
            phase,
            spin_frequency,
            ephem_vectors,
            level,
            tzr,
        )
    end
end

TOA(value, error, observing_frequency, phase, ephem_vectors) =
    TOA(value, error, observing_frequency, phase, frequency(-1.0), ephem_vectors, 0, false)

copy(toa::TOA) = TOA(
    toa.value,
    toa.error,
    toa.observing_frequency,
    toa.phase,
    toa.spin_frequency,
    toa.ephem_vectors,
    toa.level,
    toa.tzr,
)

# is_barycentered(toa::TOA) = isnothing(toa.ephem_vectors)

make_tzr_toa(tzrtdb, tzrfreq, ephem_vectors) = TOA(
    tzrtdb,
    time(0.0),
    tzrfreq,
    dimensionless(Float128(0.0)),
    frequency(-1.0),
    ephem_vectors,
    0,
    true,
)

function correct_toa_delay!(toa::TOA, delay::GQ)
    toa.value -= delay
    toa.level += 1
end

function correct_toa_phase!(toa::TOA, phase::GQ)
    toa.phase += phase
    toa.level += 1
end
