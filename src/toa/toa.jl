import Base.copy, Base.show

export TOA,
    make_tzr_toa,
    CorrectedTOA,
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
    barycentered::Bool
    tzr::Bool
    ephem::SolarSystemEphemeris
    index::UInt

    function TOA(
        value,
        error,
        observing_frequency,
        pulse_number,
        barycentered,
        tzr,
        ephem,
        index,
    )
        @assert tzr || index > 0 "The index can be 0 only for the TZR TOA."

        return new(
            value,
            error,
            observing_frequency,
            pulse_number,
            barycentered,
            tzr,
            ephem,
            index,
        )
    end
end

TOA(value, error, observing_frequency, pulse_number, barycentered, ephem, index) =
    TOA(value, error, observing_frequency, pulse_number, barycentered, false, ephem, index)

"""Create a TZR TOA object."""
make_tzr_toa(tzrtdb, tzrfreq, tzrbary, tzrephem) = TOA(
    tzrtdb,
    time(0.0),
    tzrfreq,
    dimensionless(Double64(0.0)),
    tzrbary,
    true,
    tzrephem,
    0,
)

abstract type CorrectedTOABase end

"""The accumulated timing & noise model corrections applied to a narrowband TOA."""
struct CorrectedTOA <: CorrectedTOABase
    toa::TOA
    delay::GQ{1,Float64}
    phase::GQ{0,Double64}
    efac::GQ{0,Float64}
    equad2::GQ{2,Float64}
    spin_frequency::GQ{-1,Float64}
    doppler::GQ{0,Float64}
    barycentered::Bool
    ssb_psr_pos::NTuple{3,GQ{0,Float64}}
    level::UInt

    function CorrectedTOA(
        toa,
        delay,
        phase,
        efac,
        equad2,
        spin_frequency,
        doppler,
        barycentered,
        ssb_psr_pos,
        level,
    )
        @assert abs(doppler) < dimensionless(1.0) "|doppler| must be less than 1."
        @assert spin_frequency >= frequency(0.0) "spin_frequency must either be a positive value or a default value of 0.0."

        @assert all(iszero.(ssb_psr_pos)) ||
                dot(ssb_psr_pos, ssb_psr_pos) ≈ dimensionless(1.0) "ssb_psr_pos must be a zero vector (representing pending computation) or a unit vector."

        return new(
            toa,
            delay,
            phase,
            efac,
            equad2,
            spin_frequency,
            doppler,
            barycentered,
            ssb_psr_pos,
            level,
        )
    end
end

CorrectedTOA(toa) = CorrectedTOA(
    toa,
    time(0.0),
    dimensionless(Double64(0.0)),
    dimensionless(1.0),
    GQ{2}(0.0),
    frequency(0.0),
    dimensionless(0.0),
    toa.barycentered,
    dimensionless.((0.0, 0.0, 0.0)),
    0,
)

"""Squared TOA uncertainty after applying EFAC and EQUAD."""
scaled_toa_error_sqr(ctoa::CorrectedTOA) =
    (ctoa.toa.error * ctoa.toa.error + ctoa.equad2) * ctoa.efac * ctoa.efac

"""Spin frequency in topocentric or barycentric frame, depending on the correction level.
The spin_frequency is originally in the pulsar frame."""
function doppler_shifted_spin_frequency(ctoa::CorrectedTOA)::GQ
    @assert !iszero(ctoa.spin_frequency) "The spin_frequency has not been set."
    return ctoa.spin_frequency * (1 + ctoa.doppler)
end

"""Observing frequency in the barycentric or pulsar frame, depending on the correction level.
The observing_frequency is originally in the topocentric frame."""
doppler_corrected_observing_frequency(ctoa::CorrectedTOA)::GQ =
    ctoa.toa.observing_frequency * (1 - ctoa.doppler)

"""TOA value after delay correction with 128-bit precision."""
corrected_toa_value_F128(ctoa::CorrectedTOA) = ctoa.toa.value - ctoa.delay

"""TOA value after delay correction with 64-bit precision."""
corrected_toa_value(ctoa::CorrectedTOA) = GQ{Float64}(ctoa.toa.value) - ctoa.delay

"""TOA phase residual"""
phase_residual(ctoa::CorrectedTOA) = ctoa.phase - ctoa.toa.pulse_number

"""Apply a correction to a CorrectedTOA object."""
correct_toa(
    ctoa::CorrectedTOA;
    delay::GQ = time(0.0),
    phase::GQ = dimensionless(0.0),
    efac::GQ = dimensionless(1.0),
    equad2::GQ = GQ{2}(0.0),
    delta_spin_frequency::GQ = frequency(0.0),
    doppler::GQ = dimensionless(0.0),
    barycentered::Bool = false,
    ssb_psr_pos::Union{Nothing,NTuple{3,GQ{0,Float64}}} = nothing,
) = CorrectedTOA(
    ctoa.toa,
    ctoa.delay + delay,
    ctoa.phase + phase,
    ctoa.efac * efac,
    ctoa.equad2 + equad2,
    ctoa.spin_frequency + delta_spin_frequency,
    ctoa.doppler + doppler,
    ctoa.barycentered || barycentered,
    isnothing(ssb_psr_pos) ? ctoa.ssb_psr_pos : ssb_psr_pos,
    ctoa.level + 1,
)

const day_to_s = 86400
show(io::IO, toa::TOA) = print(
    io,
    "$(toa.tzr ? "TZR" : "")TOA[MJD:$(trunc(Int, toa.value.x/day_to_s)), Freq(MHz):$(trunc(Int, toa.observing_frequency.x/1e6))]",
)
show(io::IO, ::MIME"text/plain", toa::TOABase) = show(io, toa)
show(io::IO, toas::Vector{TOA}) = print(io, "[Vector containing $(length(toas)) TOAs.]")
show(io::IO, ::MIME"text/plain", toas::Vector{T}) where {T<:TOABase} = show(io, toas)
