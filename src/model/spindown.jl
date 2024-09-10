export Spindown, phase, spin_frequency, read_params_from_dict

"""Rotation or the pulsar represented as a Taylor series in spin frequency.

Corresponds to `Spindown` in `PINT`.

Reference:
    [Backer & Hellings 1986](http://doi.org/10.1146/annurev.aa.24.090186.002541)
"""
struct Spindown <: PhaseComponent end

"""Rotational phase of the pulsar."""
function phase(::Spindown, ctoa::CorrectedTOA, params::NamedTuple)::GQ{0,Double64}
    t0 = params.PEPOCH
    t = corrected_toa_value_F128(ctoa)
    fs = params.F
    phase0 = dimensionless(0.0)
    t_t0 = t - t0
    return params.F_ * t_t0 + taylor_horner_integral(t_t0, fs, phase0)
end

"""Instantaneous rotational frequency of the pulsar."""
function spin_frequency(::Spindown, ctoa::CorrectedTOA, params::NamedTuple)::GQ{-1,Float64}
    t0 = params.PEPOCH
    t = corrected_toa_value(ctoa)
    fs = params.F
    return params.F_ + taylor_horner(t - t0, fs)
end

"""Update the `CorrectedTOA` object with the rotational phase and spin frequency"""
correct_toa(spindown::Spindown, ctoa::CorrectedTOA, params::NamedTuple) = correct_toa(
    ctoa;
    phase = phase(spindown, ctoa, params),
    delta_spin_frequency = spin_frequency(spindown, ctoa, params),
)
