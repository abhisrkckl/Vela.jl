export PowerlawRedNoiseGP,
    PowerlawDispersionNoiseGP, PowerlawChromaticNoiseGP, calc_noise_weights_inv

"""
    powerlaw(A, γ, f, f1)

Powerlaw spectral model for pulsar red noise.

Reference:
    [Lentati+ 2014](https://doi.org/10.1093/mnras/stt2122)
"""
function powerlaw(A, γ, f, f1)
    fyr = frequency(1 / 3600 / 24 / 365.25)
    denom = 12 * π^2 * fyr * fyr * fyr
    return A * A / denom * (fyr / f)^γ * f1
end

"""
    evaluate_powerlaw_red_noise_gp

Evaluate the power law Fourier-basis red noise delay/DM/CM.
"""
function evaluate_powerlaw_red_noise_gp(log10_A, γ, αs, βs, f1, Δt, ln_js)
    @assert length(αs) == length(βs) == length(ln_js)

    A = exp10(log10_A)

    ϕ1 = 2π * f1 * Δt
    exp_im_ϕ1 = exp(im * value(ϕ1))

    σ1 = sqrt(powerlaw(A, γ, f1, f1))

    result = dimensionless(0.0)

    # Handle log-spaced harmonics
    nlog = findfirst(iszero, ln_js) - 1
    for ii = 1:nlog
        α, β, ln_j = αs[ii], βs[ii], ln_js[ii]
        j = exp(ln_j)
        jfac = exp(-(γ / 2) * ln_j)
        sincosϕ = sincos(j * ϕ1)
        result += jfac * dot((α, β), sincosϕ)
    end

    # Handle linearly spaced harmonics
    ntot = length(ln_js)
    exp_im_ϕj = exp_im_ϕ1
    for ii = (nlog+1):ntot
        α, β, ln_j = αs[ii], βs[ii], ln_js[ii]
        jfac = exp(-(γ / 2) * ln_j)
        sincosϕ = imag(exp_im_ϕj), real(exp_im_ϕj)
        result += jfac * dot((α, β), sincosϕ)
        exp_im_ϕj *= exp_im_ϕ1
    end

    return σ1 * result
end

"""
    evaluate_powerlaw_red_noise_weights_inv(log10_A, γ, f1, ln_js)

Evaluate the prior weights for a power law red noise.
"""
function evaluate_powerlaw_red_noise_weights_inv(log10_A, γ, f1, ln_js)
    A = exp10(log10_A)

    P1 = powerlaw(A, γ, f1, f1)
    @assert udim(P1) == 2
    P1 = value(P1)

    npar = length(ln_js)
    weights_inv = Vector{Float64}(undef, 2 * npar)
    @inbounds for p = 1:npar
        weights_inv[p] = weights_inv[p+npar] = 1 / (P1 * exp(-γ * ln_js[p]))
    end

    return weights_inv
end

"""Compute log(j) for a given number of linear and log-spaced harmonics."""
function _calc_ln_js(Nlin, Nlog, logfac)
    log_js_lin = map(log, 1:Nlin)
    log_js_log = (1 / logfac) .^ (Nlog:-1:1)
    return vcat(log_js_log, log_js_lin)
end

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

is_gp_noise(::PowerlawRedNoiseGP) = true # coverage: ignore
get_gp_npars(arn::PowerlawRedNoiseGP) = 2 * length(arn.ln_js)

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

is_gp_noise(::PowerlawDispersionNoiseGP) = true # coverage: ignore
get_gp_npars(dmn::PowerlawDispersionNoiseGP) = 2 * length(dmn.ln_js)

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

"""
    PowerlawDispersionNoiseGP

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

is_gp_noise(::PowerlawChromaticNoiseGP) = true # coverage: ignore
get_gp_npars(crn::PowerlawChromaticNoiseGP) = 2 * length(crn.ln_js)

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
