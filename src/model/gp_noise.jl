export PowerlawRedNoiseGP, PowerlawDispersionNoiseGP, PowerlawChromaticNoiseGP

function powerlaw(A, γ, f, f1)
    fyr = frequency(1 / 3600 / 24 / 365.25)
    denom = 12 * π^2 * fyr * fyr * fyr
    return A * A / denom * (fyr / f)^γ * f1
end

function evaluate_powerlaw_red_noise_gp(log10_A, γ, αs, βs, f1, Δt, ln_js)
    @assert length(αs) == length(βs) == length(ln_js)

    A = exp10(log10_A)

    ϕ1 = 2π * f1 * Δt
    exp_im_ϕ1 = exp(im * value(ϕ1))
    exp_im_ϕjs = cumprod(map(x -> exp_im_ϕ1, ln_js))

    σ1 = sqrt(powerlaw(A, γ, f1, f1))

    result = dimensionless(0.0)
    for (α, β, ln_j, exp_im_ϕj) in zip(αs, βs, ln_js, exp_im_ϕjs)
        jfac = exp(-(γ / 2) * ln_j)
        sincosϕ = imag(exp_im_ϕj), real(exp_im_ϕj)
        result += jfac * dot((α, β), sincosϕ)
    end

    return σ1 * result
end

"""
    PowerlawRedNoiseGP

A Fourier series Gaussian process representation of the achromatic red noise where the 
power spectral density is assumed to be a power law.
"""
struct PowerlawRedNoiseGP{N} <: RedNoiseBase
    ln_js::NTuple{N,Float32}

    PowerlawRedNoiseGP(N::Int) = new{N}(Tuple(map(Float32 ∘ log, 1:N)))
end

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

"""
    PowerlawDispersionNoiseGP

A Fourier series Gaussian process representation of the dispersion noise where the 
power spectral density is assumed to be a power law.
"""
struct PowerlawDispersionNoiseGP{N} <: DispersionNoiseBase
    ln_js::NTuple{N,Float64}

    PowerlawDispersionNoiseGP(N::Int) = new{N}(Tuple(map(log, 1:N)))
end

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

struct PowerlawChromaticNoiseGP{N} <: ChromaticComponent
    ln_js::NTuple{N,Float64}

    PowerlawChromaticNoiseGP(N::Int) = new{N}(Tuple(map(log, 1:N)))
end

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
