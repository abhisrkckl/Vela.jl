export WaveX, DMWaveX

struct WaveX <: DelayComponent end

fourier_term(Φ, a, b) = dot((a, b), sincos(Φ))

function evaluate_xwavex(ctoa::CorrectedTOA, epoch::GQ, as::NTuple, bs::NTuple, fs::NTuple)
    t_t0 = corrected_toa_value(ctoa) - epoch
    k = 2 * π * t_t0
    Φs = map(f -> k * f, fs)
    return mapreduce(fourier_term, +, Φs, as, bs)
end

delay(::WaveX, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    evaluate_xwavex(ctoa, params.WXEPOCH, params.WXSIN_, params.WXCOS_, params.WXFREQ_)

struct DMWaveX <: DispersionComponent end

dispersion_slope(::DMWaveX, ctoa::CorrectedTOA, params::NamedTuple)::GQ = evaluate_xwavex(
    ctoa,
    params.DMWXEPOCH,
    params.DMWXSIN_,
    params.DMWXCOS_,
    params.DMWXFREQ_,
)

struct CMWaveX <: ChromaticComponent end

chromatic_slope(::CMWaveX, ctoa::CorrectedTOA, params::NamedTuple)::GQ = evaluate_xwavex(
    ctoa,
    params.CMWXEPOCH,
    params.CMWXSIN_,
    params.CMWXCOS_,
    params.CMWXFREQ_,
)
