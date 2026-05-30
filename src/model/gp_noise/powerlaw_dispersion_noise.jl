export PowerlawDispersionNoiseGP

"""
    PowerlawDispersionNoiseGP

A Fourier series Gaussian process representation of the dispersion noise where the 
power spectral density is assumed to be a power law. Corresponds to `PLDMNoise` in
`PINT`.

Reference:
    [Lentati+ 2014](https://doi.org/10.1093/mnras/stt2122),
    [van Haasteren & Vallisneri 2014](https://doi.org/10.1093/mnras/stu2157)
"""
struct PowerlawDispersionNoiseGP{N} <: DispersionNoiseBase
    ln_js::NTuple{N,Float64}

    PowerlawDispersionNoiseGP(Nlin::Int, Nlog::Int, logfac) =
        new{Nlin+Nlog}(Tuple(_calc_ln_js(Nlin, Nlog, logfac)))
end

is_gp_noise(::PowerlawDispersionNoiseGP) = true # COV_EXCL_LINE
get_gp_npars(dmn::PowerlawDispersionNoiseGP) = 2 * length(dmn.ln_js)
get_marginalized_param_names(dmn::PowerlawDispersionNoiseGP) =
    marginalized_param_names_for_gp_noise("PLDM", length(dmn.ln_js))

function dispersion_slope(
    dmn::PowerlawDispersionNoiseGP,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)
    νref = frequency(1.4e9)
    return (νref * νref) * evaluate_powerlaw_red_noise_gp(
        params.TNDMAMP,
        params.TNDMGAM,
        params.PLDMSIN_,
        params.PLDMCOS_,
        params.PLDMFREQ,
        corrected_toa_value(toa, toacorr, Float64) - params.PLDMEPOCH,
        dmn.ln_js,
    )
end

calc_noise_weights_inv(dmn::PowerlawDispersionNoiseGP, params::NamedTuple) =
    evaluate_powerlaw_red_noise_weights_inv(
        params.TNDMAMP,
        params.TNDMGAM,
        params.PLDMFREQ,
        dmn.ln_js,
    )

show(io::IO, dmn::PowerlawDispersionNoiseGP) =
    print(io, "PowerlawDispersionNoiseGP($(length(dmn.ln_js)) harmonics)")
