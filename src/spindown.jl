using GeometricUnits
using Quadmath

export Spindown, phase, spin_frequency, correct_toa

struct Spindown <: PhaseComponent end

_to_f128(q::GQ) = GQ(Float128(q.x), q.d)

function phase(::Spindown, toa::TOA, params)::GQ{Float128}
    t0 = params["PEPOCH"][1]
    t = toa.value
    fs = params["F"]
    phase0 = dimensionless(0.0)
    return taylor_horner_integral(t - t0, fs, phase0)
end

function spin_frequency(::Spindown, toa::TOA, params)
    t0 = params["PEPOCH"][1]
    t = toa.value
    fs = params["F"]
    f = taylor_horner(t - t0, fs)
    return quantity_like(fs[1], f.x)
end

correct_toa(spindown::Spindown, toa::TOA, params) = TOA(
    toa.value,
    toa.error,
    toa.observing_frequency,
    toa.phase + phase(spindown, toa, params),
    spin_frequency(spindown, toa, params),
    toa.ephem_vectors,
    toa.level + 1,
    toa.tzr,
)
