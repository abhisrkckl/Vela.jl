export Spindown, phase, spin_frequency, read_params_from_dict

struct Spindown <: PhaseComponent end

function phase(::Spindown, ctoa::CorrectedTOA, params::NamedTuple)::GQ{Double64}
    t0 = params.PEPOCH
    t = corrected_toa_value_F128(ctoa)
    fs = params.F
    phase0 = dimensionless(0.0)
    t_t0 = t - t0
    return params.F_ * t_t0 + taylor_horner_integral(t_t0, fs, phase0)
end

function spin_frequency(::Spindown, ctoa::CorrectedTOA, params::NamedTuple)::GQ{Float64}
    t0 = params.PEPOCH
    t = corrected_toa_value(ctoa)
    fs = params.F
    return params.F_ + taylor_horner(t - t0, fs)
end

correct_toa(spindown::Spindown, ctoa::CorrectedTOA, params::NamedTuple) = correct_toa(
    ctoa;
    phase = phase(spindown, ctoa, params),
    delta_spin_frequency = spin_frequency(spindown, ctoa, params),
)
