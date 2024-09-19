import Base.copy, Base.show

export TOA,
    make_tzr_toa,
    is_tzr,
    is_barycentered,
    TOACorrection,
    scaled_toa_error_sqr,
    doppler_shifted_spin_frequency,
    doppler_corrected_observing_frequency,
    corrected_toa_value,
    phase_residual,
    correct_toa


"""Abstract base type of all TOAs."""
abstract type TOABase end

"""A single narrowband TOA observation.

`value` incorporates the clock corrections and `ephem` contains the 
solar system ephemerides. These are computed using `PINT`.

References:
    [Hobbs+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10302.x)
    [Luo+ 2021](http://doi.org/10.3847/1538-4357/abe62f)
"""
struct TOA <: TOABase
    value::GQ{1,Double64}
    error::GQ{1,Float64}
    observing_frequency::GQ{-1,Float64}
    pulse_number::GQ{0,Double64}
    ephem::SolarSystemEphemeris
    index::UInt

    function TOA(value, error, observing_frequency, pulse_number, ephem, index)
        return new(value, error, observing_frequency, pulse_number, ephem, index)
    end
end

is_tzr(toa::TOA) = iszero(toa.index)
is_barycentered(toa::TOA) = all(iszero, toa.ephem.ssb_obs_pos)

"""Create a TZR TOA object."""
make_tzr_toa(tzrtdb, tzrfreq, tzrephem) =
    TOA(tzrtdb, time(0.0), tzrfreq, dimensionless(Double64(0.0)), tzrephem, 0)

abstract type TOACorrectionBase end

"""The accumulated timing & noise model corrections applied to a narrowband TOA."""
struct TOACorrection <: TOACorrectionBase
    delay::GQ{1,Float64}
    phase::GQ{0,Double64}
    efac::GQ{0,Float64}
    equad2::GQ{2,Float64}
    spin_frequency::GQ{-1,Float64}
    doppler::GQ{0,Float64}
    ssb_psr_pos::NTuple{3,GQ{0,Float64}}

    function TOACorrection(delay, phase, efac, equad2, spin_frequency, doppler, ssb_psr_pos)
        @assert abs(doppler) < dimensionless(1.0) "|doppler| must be less than 1."
        @assert spin_frequency >= frequency(0.0) "spin_frequency must either be a positive value or a default value of 0.0."

        @assert all(iszero.(ssb_psr_pos)) ||
                dot(ssb_psr_pos, ssb_psr_pos) ≈ dimensionless(1.0) "ssb_psr_pos must be a zero vector (representing pending computation) or a unit vector."

        return new(delay, phase, efac, equad2, spin_frequency, doppler, ssb_psr_pos)
    end
end

is_barycentered(toacorr::TOACorrection) = !all(iszero, toacorr.ssb_psr_pos)

TOACorrection() = TOACorrection(
    time(0.0),
    dimensionless(Double64(0.0)),
    dimensionless(1.0),
    GQ{2}(0.0),
    frequency(0.0),
    dimensionless(0.0),
    dimensionless.((0.0, 0.0, 0.0)),
)

"""Squared TOA uncertainty after applying EFAC and EQUAD."""
scaled_toa_error_sqr(toa::TOA, toacorr::TOACorrection) =
    (toa.error * toa.error + toacorr.equad2) * toacorr.efac * toacorr.efac

"""Spin frequency in topocentric or barycentric frame, depending on the correction level.
The spin_frequency is originally in the pulsar frame."""
function doppler_shifted_spin_frequency(toacorr::TOACorrection)
    @assert !iszero(toacorr.spin_frequency) "The spin_frequency has not been set."
    return toacorr.spin_frequency * (1 + toacorr.doppler)
end

"""Observing frequency in the barycentric or pulsar frame, depending on the correction level.
The observing_frequency is originally in the topocentric frame."""
doppler_corrected_observing_frequency(toa::TOA, toacorr::TOACorrection) =
    toa.observing_frequency * (1 - toacorr.doppler)

"""TOA value after delay correction."""
corrected_toa_value(toa::TOA, toacorr::TOACorrection)::GQ{1,Double64} =
    toa.value - toacorr.delay
corrected_toa_value(toa::TOA, toacorr::TOACorrection, ::Type{Float64}) =
    GQ{Float64}(corrected_toa_value(toa, toacorr))

"""TOA phase residual"""
phase_residual(toa::TOA, toacorr::TOACorrection) = toacorr.phase - toa.pulse_number

"""Apply a correction to a CorrectedTOA object."""
correct_toa(
    toacorr::TOACorrection;
    delay::GQ = time(0.0),
    phase::GQ = dimensionless(0.0),
    efac::GQ = dimensionless(1.0),
    equad2::GQ = GQ{2}(0.0),
    delta_spin_frequency::GQ = frequency(0.0),
    doppler::GQ = dimensionless(0.0),
    ssb_psr_pos::Union{Nothing,NTuple{3,GQ{0,Float64}}} = nothing,
) = TOACorrection(
    toacorr.delay + delay,
    toacorr.phase + phase,
    toacorr.efac * efac,
    toacorr.equad2 + equad2,
    toacorr.spin_frequency + delta_spin_frequency,
    toacorr.doppler + doppler,
    isnothing(ssb_psr_pos) ? toacorr.ssb_psr_pos : ssb_psr_pos,
)

const day_to_s = 86400
show(io::IO, toa::TOA) = print(
    io,
    "TOA[MJD:$(trunc(Int, toa.value.x/day_to_s)), Freq(MHz):$(trunc(Int, toa.observing_frequency.x/1e6))]",
)
show(io::IO, ::MIME"text/plain", toa::TOABase) = show(io, toa)
show(io::IO, toas::Vector{TOA}) = print(io, "[Vector containing $(length(toas)) TOAs.]")
show(io::IO, ::MIME"text/plain", toas::Vector{T}) where {T<:TOABase} = show(io, toas)
