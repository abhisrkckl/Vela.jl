export WaveX, DMWaveX, CMWaveX

"""A Fourier series representation of the achromatic red noise.

Reference:
    [Susobhanan+ 2024](http://doi.org/10.3847/1538-4357/ad59f7)
"""
struct WaveX <: DelayComponent end

"""A single term in the Fourier series appearing in the implementation of `WaveX`,
`DMWaveX`, and `CMWaveX`."""
fourier_term(Φ, a, b) = dot((a, b), sincos(Φ))

"""Evaluate the Fourier series appearing in the implementation of `WaveX`,
`DMWaveX`, and `CMWaveX`."""
function evaluate_xwavex(
    toa::TOA,
    toacorr::TOACorrection,
    epoch::GQ,
    as::NTuple,
    bs::NTuple,
    fs::NTuple,
)
    t_t0 = corrected_toa_value(toa, toacorr, Float64) - epoch
    k = 2 * π * t_t0
    Φs = map(f -> k * f, fs)
    return mapreduce(fourier_term, +, Φs, as, bs)
end

"""Delay due to achromatic red noise (Fourier series representation)."""
delay(::WaveX, toa::TOA, toacorr::TOACorrection, params::NamedTuple)::GQ = evaluate_xwavex(
    toa,
    toacorr,
    params.WXEPOCH,
    params.WXSIN_,
    params.WXCOS_,
    params.WXFREQ_,
)

"""A Fourier series representation of the achromatic red noise.

Reference:
    [Susobhanan+ 2024](http://doi.org/10.3847/1538-4357/ad59f7)
"""
struct DMWaveX <: DispersionComponent end

dispersion_slope(::DMWaveX, toa::TOA, toacorr::TOACorrection, params::NamedTuple)::GQ =
    evaluate_xwavex(
        toa,
        toacorr,
        params.DMWXEPOCH,
        params.DMWXSIN_,
        params.DMWXCOS_,
        params.DMWXFREQ_,
    )

"""A Fourier series representation of the variable-index chromatic red noise.

Reference:
    [Susobhanan+ 2024](http://doi.org/10.3847/1538-4357/ad59f7)
"""
struct CMWaveX <: ChromaticComponent end

"""Chromatic slope due to variable-index chromatic red noise (Fourier series representation)."""
chromatic_slope(::CMWaveX, toa::TOA, toacorr::TOACorrection, params::NamedTuple)::GQ =
    evaluate_xwavex(
        toa,
        toacorr,
        params.CMWXEPOCH,
        params.CMWXSIN_,
        params.CMWXCOS_,
        params.CMWXFREQ_,
    )
