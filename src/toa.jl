using GeometricUnits
using Quadmath
using LinearAlgebra
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

"""A single narrowband TOA observation."""
struct TOA
    value::GQ{Float128}
    error::GQ{Float64}
    observing_frequency::GQ{Float64}
    pulse_number::GQ{Float128}
    barycentered::Bool
    tzr::Bool
    ephem::SolarSystemEphemeris

    function TOA(value, error, observing_frequency, pulse_number, barycentered, tzr, ephem)
        @assert value.d == 1 "Dimension mismatch in value (given $(value.d), expected 1)."
        @assert error.d == 1 "Dimension mismatch in error (given $(error.d), expected 1)."
        @assert observing_frequency.d == -1 "Dimension mismatch in observing_frequency (given $(observing_frequency.d), expected -1)."
        @assert pulse_number.d == 0 "Dimension mismatch in pulse_number (given $(pulse_number.d), expected 0)."

        return new(
            value,
            error,
            observing_frequency,
            pulse_number,
            barycentered,
            tzr,
            ephem,
        )
    end
end

TOA(value, error, observing_frequency, pulse_number, barycentered, ephem) =
    TOA(value, error, observing_frequency, pulse_number, barycentered, false, ephem)

"""Create a TZR TOA object."""
make_tzr_toa(tzrtdb, tzrfreq, tzrbary, tzrephem) =
    TOA(tzrtdb, time(0.0), tzrfreq, dimensionless(Float128(0.0)), tzrbary, true, tzrephem)

"""The accumulated timing & noise model corrections applied to a narrowband TOA."""
struct CorrectedTOA
    toa::TOA
    delay::GQ{Float64}
    phase::GQ{Float128}
    efac::GQ{Float64}
    equad2::GQ{Float64}
    spin_frequency::GQ{Float64}
    doppler::GQ{Float64}
    barycentered::Bool
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
        level,
    )
        @assert delay.d == 1 "Dimension mismatch in value (given $(delay.d), expected 1)."
        @assert phase.d == 0 "Dimension mismatch in phase (given $(phase.d), expected 0)."
        @assert efac.d == 0 "Dimension mismatch in efac (given $(efac.d), expected 0)."
        @assert equad2.d == 2 "Dimension mismatch in equad (given $(equad2.d), expected 2)."
        @assert spin_frequency.d == -1 "Dimension mismatch in spin_frequency (given $(spin_frequency.d), expected -1)."
        @assert doppler.d == 0 "Dimension mismatch in doppler (given $(doppler.d), expected 0)."

        @assert abs(doppler) < dimensionless(1.0) "|doppler| must be less than 1."
        @assert spin_frequency == frequency(-1.0) || spin_frequency > frequency(0.0) "spin_frequency must either be a positive value or a default value of -1."

        return new(
            toa,
            delay,
            phase,
            efac,
            equad2,
            spin_frequency,
            doppler,
            barycentered,
            level,
        )
    end
end

CorrectedTOA(toa) = CorrectedTOA(
    toa,
    time(0.0),
    dimensionless(Float128(0.0)),
    dimensionless(1.0),
    GQ(0.0, 2),
    frequency(-1.0),
    dimensionless(0.0),
    toa.barycentered,
    0,
)

"""Squared TOA uncertainty after applying EFAC and EQUAD."""
scaled_toa_error_sqr(ctoa::CorrectedTOA) =
    (ctoa.toa.error * ctoa.toa.error + ctoa.equad2) * ctoa.efac * ctoa.efac

"""Spin frequency in topocentric or barycentric frame, depending on the correction level.
The spin_frequency is originally in the pulsar frame."""
function doppler_shifted_spin_frequency(ctoa::CorrectedTOA)::GQ
    @assert ctoa.spin_frequency != frequency(-1.0) "The spin_frequency has not been set."
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
    equad2::GQ = GQ(0.0, 2),
    delta_spin_frequency::GQ = frequency(0.0),
    doppler::GQ = dimensionless(0.0),
    barycentered = false,
) = CorrectedTOA(
    ctoa.toa,
    ctoa.delay + delay,
    ctoa.phase + phase,
    ctoa.efac * efac,
    ctoa.equad2 + equad2,
    ctoa.spin_frequency + delta_spin_frequency,
    ctoa.doppler + doppler,
    ctoa.barycentered || barycentered,
    ctoa.level + 1,
)

const day_to_s = 86400
show(io::IO, toa::TOA) = print(
    io,
    "$(toa.tzr ? "TZR" : "")TOA[MJD:$(trunc(Int, toa.value.x/day_to_s)), Freq(MHz):$(trunc(Int, toa.observing_frequency.x/1e6))]",
)
show(io::IO, ::MIME"text/plain", toa::TOA) = show(io, toa)
show(io::IO, toas::Vector{TOA}) = print(io, "[Vector containing $(length(toas)) TOAs.]")
show(io::IO, ::MIME"text/plain", toas::Vector{TOA}) = show(io, toas)
