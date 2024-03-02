using GeometricUnits
using Quadmath

export Spindown, phase, spin_frequency, read_params_from_dict

struct Spindown <: PhaseComponent end

read_params_from_dict(::Spindown, params::Dict) =
    (PEPOCH = params["PEPOCH"][1], F = params["F"])

function phase(::Spindown, toa::TOA, params::NamedTuple)::GQ{Float128}
    t0 = params.PEPOCH
    t = toa.value
    fs = params.F
    phase0 = dimensionless(0.0)
    return @fastmath(taylor_horner_integral(t - t0, fs, phase0))
end

function spin_frequency(::Spindown, toa::TOA, params::NamedTuple)
    t0 = params.PEPOCH
    t = GQ{Float64}(toa.value)
    fs = params.F
    return taylor_horner(t - t0, fs)
end

correct_toa(spindown::Spindown, toa::TOA, params::NamedTuple) = TOA(
    toa.value,
    toa.error,
    toa.observing_frequency,
    toa.phase + phase(spindown, toa, params),
    spin_frequency(spindown, toa, params),
    toa.barycentered,
    toa.tzr,
    toa.level + 1,
)
