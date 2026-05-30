export PowerlawChromaticNoiseGP

"""
    PowerlawChromaticNoiseGP

A Fourier series Gaussian process representation of the chromatic noise where the 
power spectral density is assumed to be a power law. Corresponds to `PLChromNoise`
in PINT.

Reference:
    [Lentati+ 2014](https://doi.org/10.1093/mnras/stt2122),
    [van Haasteren & Vallisneri 2014](https://doi.org/10.1093/mnras/stu2157)
"""
struct PowerlawChromaticNoiseGP{N} <: ChromaticNoiseBase
    ln_js::NTuple{N,Float64}

    PowerlawChromaticNoiseGP(Nlin::Int, Nlog::Int, logfac) =
        new{Nlin+Nlog}(Tuple(_calc_ln_js(Nlin, Nlog, logfac)))
end

is_gp_noise(::PowerlawChromaticNoiseGP) = true # COV_EXCL_LINE
get_gp_npars(crn::PowerlawChromaticNoiseGP) = 2 * length(crn.ln_js)
get_marginalized_param_names(crn::PowerlawChromaticNoiseGP) =
    marginalized_param_names_for_gp_noise("PLCHROM", length(crn.ln_js))

function chromatic_slope(
    crn::PowerlawChromaticNoiseGP,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)
    νref_val = 1400.0
    return (νref_val^params.TNCHROMIDX) * evaluate_powerlaw_red_noise_gp(
        params.TNCHROMAMP,
        params.TNCHROMGAM,
        params.PLCHROMSIN_,
        params.PLCHROMCOS_,
        params.PLCHROMFREQ,
        corrected_toa_value(toa, toacorr, Float64) - params.PLCHROMEPOCH,
        crn.ln_js,
    )
end

calc_noise_weights_inv(crn::PowerlawChromaticNoiseGP, params::NamedTuple) =
    evaluate_powerlaw_red_noise_weights_inv(
        params.TNCHROMAMP,
        params.TNCHROMGAM,
        params.PLCHROMFREQ,
        crn.ln_js,
    )

show(io::IO, crn::PowerlawChromaticNoiseGP) =
    print(io, "PowerlawChromaticNoiseGP($(length(crn.ln_js)) harmonics)")
