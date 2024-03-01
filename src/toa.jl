using GeometricUnits
using Quadmath
using LinearAlgebra
import Base.copy, Base.show

export TOA, correct_toa_delay, correct_toa_phase, make_tzr_toa

struct TOA
    value::GQ{Float128}
    error::GQ{Float64}
    observing_frequency::GQ{Float64}
    phase::GQ{Float128}
    spin_frequency::GQ{Float64}
    barycentered::Bool
    tzr::Bool
    level::UInt

    function TOA(
        value,
        error,
        observing_frequency,
        phase,
        spin_frequency,
        barycentered,
        tzr,
        level,
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
            barycentered,
            tzr,
            level,
        )
    end
end

TOA(value, error, observing_frequency, phase, barycentered) =
    TOA(value, error, observing_frequency, phase, frequency(-1.0), barycentered, false, 0)

copy(toa::TOA) = TOA(
    toa.value,
    toa.error,
    toa.observing_frequency,
    toa.phase,
    toa.spin_frequency,
    toa.barycentered,
    toa.tzr,
    toa.level,
)

make_tzr_toa(tzrtdb, tzrfreq, tzrbary) = TOA(
    tzrtdb,
    time(0.0),
    tzrfreq,
    dimensionless(Float128(0.0)),
    frequency(-1.0),
    tzrbary,
    true,
    0,
)

correct_toa_delay(toa::TOA, delay::GQ) = TOA(
    toa.value - delay,
    toa.error,
    toa.observing_frequency,
    toa.phase,
    toa.spin_frequency,
    toa.barycentered,
    toa.tzr,
    toa.level + 1,
)

correct_toa_phase(toa::TOA, phase::GQ) = TOA(
    toa.value,
    toa.error,
    toa.observing_frequency,
    toa.phase + phase,
    toa.spin_frequency,
    toa.barycentered,
    toa.tzr,
    toa.level + 1,
)

const day_to_s = 86400
const tzrstr = "TZR"
show(io::IO, toa::TOA) = print(
    io,
    "$(toa.tzr ? "TZR" : "")TOA[MJD:$(trunc(Int, toa.value.x/day_to_s)), Freq(MHz):$(trunc(Int, toa.observing_frequency.x/1e6))]",
)
show(io::IO, ::MIME"text/plain", toa::TOA) = show(io, toa)
show(io::IO, toas::Vector{TOA}) = print(io, "[Vector of $(length(toas)) TOAs.]")
show(io::IO, ::MIME"text/plain", toas::Vector{TOA}) = show(io, toas)
