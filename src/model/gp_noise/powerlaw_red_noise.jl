export PowerlawRedNoiseGP, calc_noise_weights_inv

"""
    PowerlawRedNoiseGP

A Fourier series Gaussian process representation of the achromatic red noise where the 
power spectral density is assumed to be a power law. Corresponds to `PLRedNoise` in 
`PINT`.

Reference:
    [Lentati+ 2014](https://doi.org/10.1093/mnras/stt2122),
    [van Haasteren & Vallisneri 2014](https://doi.org/10.1093/mnras/stu2157)
"""
struct PowerlawRedNoiseGP{N} <: RedNoiseBase
    ln_js::NTuple{N,Float32}

    PowerlawRedNoiseGP(Nlin::Int, Nlog::Int, logfac::Float64) =
        new{Nlin+Nlog}(Tuple(_calc_ln_js(Nlin, Nlog, logfac)))
end

is_gp_noise(::PowerlawRedNoiseGP) = true # COV_EXCL_LINE
get_gp_npars(arn::PowerlawRedNoiseGP) = 2 * length(arn.ln_js)
get_marginalized_param_names(arn::PowerlawRedNoiseGP) =
    marginalized_param_names_for_gp_noise("PLRED", length(arn.ln_js))

delay(arn::PowerlawRedNoiseGP, toa::TOA, toacorr::TOACorrection, params::NamedTuple) =
    evaluate_powerlaw_red_noise_gp(
        params.TNREDAMP,
        params.TNREDGAM,
        params.PLREDSIN_,
        params.PLREDCOS_,
        params.PLREDFREQ,
        corrected_toa_value(toa, toacorr, Float64) - params.PLREDEPOCH,
        arn.ln_js,
    )

calc_noise_weights_inv(arn::PowerlawRedNoiseGP, params::NamedTuple) =
    evaluate_powerlaw_red_noise_weights_inv(
        params.TNREDAMP,
        params.TNREDGAM,
        params.PLREDFREQ,
        arn.ln_js,
    )

show(io::IO, arn::PowerlawRedNoiseGP) =
    print(io, "PowerlawRedNoiseGP($(length(arn.ln_js)) harmonics)")
