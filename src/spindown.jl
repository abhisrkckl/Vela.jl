using GeometricUnits
using Quadmath

export Spindown, phase, spin_frequency, read_params_from_dict

struct Spindown <: PhaseComponent end

function phase(::Spindown, ctoa::CorrectedTOA, params::NamedTuple)::GQ{Float128}
    t0 = params.PEPOCH
    t = corrected_toa_value_F128(ctoa)
    fs = params.F
    phase0 = dimensionless(0.0)
    return taylor_horner_integral(t - t0, fs, phase0)
end

function spin_frequency(::Spindown, ctoa::CorrectedTOA, params::NamedTuple)::GQ{Float64}
    t0 = params.PEPOCH
    t = corrected_toa_value(ctoa)
    fs = params.F
    return taylor_horner(t - t0, fs)
end

correct_toa(spindown::Spindown, ctoa::CorrectedTOA, params::NamedTuple) = correct_toa(
    ctoa;
    phase = phase(spindown, ctoa, params),
    delta_spin_frequency = spin_frequency(spindown, ctoa, params),
)

const priors_F = (
    LogUniform(1e-4, 1e+3),
    Uniform(-1e-9, 1e-9),
    Uniform(-1e-18, 1e-18),
    Uniform(-1e-27, 1e-27),
    Uniform(-1e-35, 1e-35),
)

"""
Prior distributions for F0 ... F4 based on psrcat values, but broadened to 
include one extra order of magnitude in either direction (except F0, which has a 
theoretical upper limit).
"""
lnprior(::Spindown, params::NamedTuple) = _lnprior(priors_F, params, :F)
