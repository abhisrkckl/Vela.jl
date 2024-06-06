using GeometricUnits
using Quadmath
using LinearAlgebra
import Base.copy, Base.show

export TOA, correct_toa_delay, correct_toa_phase, make_tzr_toa

"""
    A type representing a single narrow-band TOA.
"""
struct TOA
    value::GQ{Float128}
    error::GQ{Float64}
    observing_frequency::GQ{Float64}
    phase::GQ{Float128}
    spin_frequency::GQ{Float64}
    doppler::GQ{Float64}
    barycentered::Bool
    tzr::Bool
    level::UInt
    ephem::SolarSystemEphemeris

    function TOA(
        value,
        error,
        observing_frequency,
        phase,
        spin_frequency,
        doppler,
        barycentered,
        tzr,
        level,
        ephem,
    )
        @assert value.d == 1 "Dimension mismatch in value (given $(value.d), expected 1)."
        @assert error.d == 1 "Dimension mismatch in error (given $(error.d), expected 1)."
        @assert observing_frequency.d == -1 "Dimension mismatch in observing_frequency (given $(observing_frequency.d), expected -1)."
        @assert spin_frequency.d == -1 "Dimension mismatch in spin_frequency (given $(spin_frequency.d), expected -1)."
        @assert phase.d == 0 "Dimension mismatch in phase (given $(phase.d), expected 0)."

        return new(
            value,
            error,
            observing_frequency,
            phase,
            spin_frequency,
            doppler,
            barycentered,
            tzr,
            level,
            ephem,
        )
    end
end

TOA(value, error, observing_frequency, phase, barycentered, ephem) = TOA(
    value,
    error,
    observing_frequency,
    phase,
    frequency(-1.0),
    dimensionless(0.0),
    barycentered,
    false,
    0,
    ephem,
)

"""Create a TZR TOA object."""
make_tzr_toa(tzrtdb, tzrfreq, tzrbary, tzrephem) = TOA(
    tzrtdb,
    time(0.0),
    tzrfreq,
    dimensionless(Float128(0.0)),
    frequency(-1.0),
    dimensionless(0.0),
    tzrbary,
    true,
    0,
    tzrephem,
)

"""Correct the TOA object by a given time delay."""
correct_toa_delay(toa::TOA, delay::GQ) = TOA(
    toa.value - delay,
    toa.error,
    toa.observing_frequency,
    toa.phase,
    toa.spin_frequency,
    toa.doppler,
    toa.barycentered,
    toa.tzr,
    toa.level + 1,
    toa.ephem,
)

"""Correct a TOA object by a given phase."""
correct_toa_phase(toa::TOA, phase::GQ) = TOA(
    toa.value,
    toa.error,
    toa.observing_frequency,
    toa.phase + phase,
    toa.spin_frequency,
    toa.doppler,
    toa.barycentered,
    toa.tzr,
    toa.level + 1,
    toa.ephem,
)

const day_to_s = 86400
show(io::IO, toa::TOA) = print(
    io,
    "$(toa.tzr ? "TZR" : "")TOA[MJD:$(trunc(Int, toa.value.x/day_to_s)), Freq(MHz):$(trunc(Int, toa.observing_frequency.x/1e6))]",
)
show(io::IO, ::MIME"text/plain", toa::TOA) = show(io, toa)
show(io::IO, toas::Vector{TOA}) = print(io, "[Vector containing $(length(toas)) TOAs.]")
show(io::IO, ::MIME"text/plain", toas::Vector{TOA}) = show(io, toas)
