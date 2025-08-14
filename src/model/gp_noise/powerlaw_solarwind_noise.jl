export PowerlawSolarWindNoiseGP

"""
    PowerlawSolarWindNoiseGP

A Fourier series Gaussian process representation of the solar wind dispersion noise 
where the power spectral density is assumed to be a power law. Corresponds to `PLSWNoise` 
in `PINT`.

Reference:
    [Susarla+ 2024](https://doi.org/10.1051/0004-6361/202450680),
"""
struct PowerlawSolarWindNoiseGP{N} <: SolarWind
    ln_js::NTuple{N,Float64}

    PowerlawSolarWindNoiseGP(Nlin::Int, Nlog::Int, logfac) =
        new{Nlin+Nlog}(Tuple(_calc_ln_js(Nlin, Nlog, logfac)))
end

is_gp_noise(::PowerlawSolarWindNoiseGP) = true # COV_EXCL_LINE
get_gp_npars(swn::PowerlawSolarWindNoiseGP) = 2 * length(swn.ln_js)
get_marginalized_param_names(swn::PowerlawSolarWindNoiseGP) =
    marginalized_param_names_for_gp_noise("PLSW", length(swn.ln_js))

function dispersion_slope(
    swn::PowerlawSolarWindNoiseGP,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)
    # This following ugly unit hacks are due to the solar wind noise
    # expression used in Susarla+ 2024 not being dimensionally correct.
    # These expressions are implemented in enterprise_extensions and PINT.

    τref = time(1.0)
    cm = distance(0.01)
    ne_sw_ref = 1/(cm*cm*cm)

    ρ, r = sun_angle_and_distance(toa, toacorr)
    t = corrected_toa_value(toa, toacorr, Float64)
    DM_sw_ref = ne_sw_ref * AU * AU * ρ / (r * sin(ρ))

    return (DM_sw_ref / τref) * evaluate_powerlaw_red_noise_gp(
        params.TNSWAMP,
        params.TNSWGAM,
        params.PLSWSIN_,
        params.PLSWCOS_,
        params.PLSWFREQ,
        t - params.PLSWEPOCH,
        swn.ln_js,
    )
end

calc_noise_weights_inv(swn::PowerlawSolarWindNoiseGP, params::NamedTuple) =
    evaluate_powerlaw_red_noise_weights_inv(
        params.TNSWAMP,
        params.TNSWGAM,
        params.PLSWFREQ,
        swn.ln_js,
    )

show(io::IO, swn::PowerlawSolarWindNoiseGP) =
    print(io, "PowerlawSolarWindNoiseGP($(length(swn.ln_js)) harmonics)")
