export Spindown, phase, spin_frequency, read_params_from_dict

struct Spindown <: PhaseComponent
    F_::GQ{Float64}
    PEPOCH::GQ{Float64}
end

_to_f128(q::GQ) = GQ(Float128(value(q)), udim(q))

function phase(sd::Spindown, ctoa::CorrectedTOA, params::NamedTuple)::GQ{Float128}
    t = corrected_toa_value_F128(ctoa)
    t0 = _to_f128(sd.PEPOCH)
    f_ = _to_f128(sd.F_)
    fs = _to_f128.(params.F)
    fs = (f_ + fs[1], fs[2:end]...)
    phase0 = dimensionless(0.0)
    t_t0 = t - t0
    return taylor_horner_integral(t_t0, fs, phase0)
end

function spin_frequency(sd::Spindown, ctoa::CorrectedTOA, params::NamedTuple)::GQ{Float64}
    t0 = sd.PEPOCH
    t = corrected_toa_value(ctoa)
    fs = params.F
    return sd.F_ + taylor_horner(t - t0, fs)
end

correct_toa(spindown::Spindown, ctoa::CorrectedTOA, params::NamedTuple) = correct_toa(
    ctoa;
    phase = phase(spindown, ctoa, params),
    delta_spin_frequency = spin_frequency(spindown, ctoa, params),
)
