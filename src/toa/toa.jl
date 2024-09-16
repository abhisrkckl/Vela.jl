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


"""
    TOABase

Abstract base type of all TOAs."""
abstract type TOABase end

"""
    TOA

A single narrowband TOA observation.

`value` is the TOA in the TDB timescale incorporating the clock corrections.
`ephem` contains the solar system ephemerides. These are computed using `PINT`.

References:
    [Hobbs+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10302.x),
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

"""
    make_tzr_toa(tzrtdb, tzrfreq, tzrbary, tzrephem)

Create a TZR TOA object.

# Arguments
- `tzrtdb::GQ{1,Double64}`: The TZR TOA value (`TZRMJD`) in the TDB timescale
- `tzrfreq::GQ{-1,Float64}`: The TZR TOA observing frequency (`TZRFRQ`)
- `tzrbary::Bool`: Whether the TZR TOA has been barycentered
- `tzrephem::SolarSystemEphemeris`: The solar system ephemeris at `TZRMJD`
"""
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

"""
    CorrectedTOABase

The abstract base type representing the accumulated timing & noise model 
corrections applied to a TOA.
"""
abstract type CorrectedTOABase end

"""
    CorrectedTOA

The accumulated timing & noise model corrections applied to a narrowband TOA."""
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

"""
    scaled_toa_error_sqr(ctoa::CorrectedTOA)

Squared TOA uncertainty incorporating the effect of `EFAC`s and `EQUAD`s.
"""
scaled_toa_error_sqr(ctoa::CorrectedTOA) =
    (ctoa.toa.error * ctoa.toa.error + ctoa.equad2) * ctoa.efac * ctoa.efac

"""
    doppler_shifted_spin_frequency(ctoa::CorrectedTOA)

The spin frequency after applying the various Doppler corrections (i.e., from the
solar system motion, pulsar binary motion, etc.) depending on the correction level. 
The spin frequency is originally computed in the pulsar frame.
"""
function doppler_shifted_spin_frequency(ctoa::CorrectedTOA)
    @assert !iszero(ctoa.spin_frequency) "The spin_frequency has not been set."
    return ctoa.spin_frequency * (1 + ctoa.doppler)
end

"""
    doppler_corrected_observing_frequency(ctoa::CorrectedTOA)

The observing radio frequency after applying the various Doppler corrections (i.e., 
from the solar system motion, pulsar binary motion, etc.) depending on the correction 
level. The observing frequency is originally measured in the topocentric frame.
"""
doppler_corrected_observing_frequency(ctoa::CorrectedTOA)::GQ =
    ctoa.toa.observing_frequency * (1 - ctoa.doppler)

"""
    corrected_toa_value_F128(ctoa::CorrectedTOA)

TOA value after delay correction with extended precision.
"""
corrected_toa_value_F128(ctoa::CorrectedTOA) = ctoa.toa.value - ctoa.delay

"""
    corrected_toa_value(ctoa::CorrectedTOA)

TOA value after delay correction with 64-bit precision.
"""
corrected_toa_value(ctoa::CorrectedTOA) = GQ{Float64}(ctoa.toa.value) - ctoa.delay

"""TOA phase residual"""
phase_residual(ctoa::CorrectedTOA) = ctoa.phase - ctoa.toa.pulse_number

"""
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
    )

Apply a correction to a TOA.

# Arguments
- `ctoa`: A `CorrectedTOA` object to be corrected
- `delay`: A delay to be added to existing delay corrections
- `phase`: A phase to be added to existing phase corrections
- `efac`: An EFAC to be multiplied to existing EFAC corrections
- `equad2`: A squared EQUAD to be added to existing squared EQUAD corrections
- `delta_spin_frequency`: A spin frequency correction (additive)
- `doppler`: A Doppler shift (additive)
- `barycentered`: Whether the TOA has been barycentered (will be OR-ed to the existing value).
- `ssb_psr_pos`: A unit 3-vector pointing from the SSB to the pulsar (will replace the existing value).
"""
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
