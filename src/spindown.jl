using GeometricUnits
using Quadmath
using .Threads

export Spindown, phase, spin_frequency, read_params_from_dict

struct Spindown <: PhaseComponent end

read_params_from_dict(::Spindown, params::Dict) =
    (PEPOCH = params["PEPOCH"][1], F = params["F"])

function phase(::Spindown, toa::TOA, params::NamedTuple)::GQ{Float128}
    t0 = params.PEPOCH
    t = toa.value
    fs = params.F
    phase0 = dimensionless(0.0)
    return taylor_horner_integral(t - t0, fs, phase0)
end

function spin_frequency(::Spindown, toa::TOA, params::NamedTuple)
    t0 = params.PEPOCH
    t = time(Float64(toa.value.x))
    fs = params.F
    f = taylor_horner(t - t0, fs)
    return quantity_like(fs[1], f.x)
end

function correct_toa!(spindown::Spindown, toa::TOA, params::NamedTuple)
    toa.phase += phase(spindown, toa, params)
    toa.spin_frequency = spin_frequency(spindown, toa, params)
    toa.level += 1
end
