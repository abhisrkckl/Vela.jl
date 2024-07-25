export WaveX, DMWaveX

struct WaveX <: DelayComponent end

fourier_term(Φ, a, b) = dot((a, b), sincos(Φ))

function delay(::WaveX, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    t_t0 = corrected_toa_value(ctoa) - params.WXEPOCH
    as, bs = params.WXSIN_, params.WXCOS_
    fs = params.WXFREQ_
    k = 2 * π * t_t0
    Φs = map(f -> k * f, fs)
    return mapreduce(fourier_term, +, Φs, as, bs)
end

struct DMWaveX <: DispersionComponent end

function dispersion_slope(::DMWaveX, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    t_t0 = corrected_toa_value(ctoa) - params.DMWXEPOCH
    as, bs = params.DMWXSIN_, params.DMWXCOS_
    fs = params.DMWXFREQ_
    k = 2 * π * t_t0
    Φs = map(f -> k * f, fs)
    return mapreduce(fourier_term, +, Φs, as, bs)
end
