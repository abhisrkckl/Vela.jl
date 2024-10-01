export PowerlawRedNoiseGP, PowerlawDispersionNoiseGP, PowerlawChromaticNoiseGP

function powerlaw(A, γ, f, f1)
    fyr = frequency(1 / 3600 / 24 / 365.25)
    denom = 12 * π^2 * fyr * fyr * fyr
    return A * A / denom * (fyr / f)^γ * f1
end

function evaluate_powerlaw_red_noise_gp(
    log10_A,
    γ,
    αs,
    βs,
    f1,
    Δt,
    ln_js,
    unit_conversion_factor = 1,
)
    @assert length(αs) == length(βs) == length(ln_js)

    Nharms = length(ln_js)

    A = 10^log10_A

    ϕ1 = 2π * f1 * Δt
    exp_im_ϕ1 = exp(im * value(ϕ1))
    exp_im_ϕjs = cumprod(map(x -> exp_im_ϕ1, ln_js))

    σ1 = sqrt(powerlaw(A, γ, f1, f1))

    result = time(0.0)
    for ii = 1:Nharms
        σ = σ1 * exp(-(γ / 2) * ln_js[ii])

        a = σ * αs[ii]
        b = σ * βs[ii]

        sincosϕ = imag(exp_im_ϕjs[ii]), real(exp_im_ϕjs[ii])
        result += dot((a, b), sincosϕ) * unit_conversion_factor
    end

    return result
end

struct PowerlawRedNoiseGP{N} <: DelayComponent
    ln_js::NTuple{N,Float32}

    PowerlawRedNoiseGP(N::Int) = new{N}(Tuple(map(Float32 ∘ log, 1:N)))
end

delay(arn::PowerlawRedNoiseGP, toa::TOA, toacorr::TOACorrection, params::NamedTuple) =
    evaluate_powerlaw_red_noise_gp(
        GQ{Float32}(params.TNREDAMP),
        GQ{Float32}(params.TNREDGAM),
        map(GQ{Float32}, params.PLREDSIN_),
        map(GQ{Float32}, params.PLREDCOS_),
        GQ{Float32}(params.PLREDFREQ),
        GQ{Float32}(corrected_toa_value(toa, toacorr, Float64) - params.PLREDEPOCH),
        arn.ln_js,
    )

struct PowerlawDispersionNoiseGP{N} <: DispersionComponent
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
    return evaluate_powerlaw_red_noise_gp(
        params.TNDMAMP,
        params.TNDMGAM,
        params.PLDMSIN_,
        params.PLDMCOS_,
        params.PLDMFREQ,
        corrected_toa_value(toa, toacorr, Float64) - params.PLDMEPOCH,
        dmn.ln_js,
        νref * νref,
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
    return evaluate_powerlaw_red_noise_gp(
        params.TNCHROMAMP,
        params.TNCHROMGAM,
        params.PLCHROMSIN_,
        params.PLCHROMCOS_,
        params.PLCHROMFREQ,
        corrected_toa_value(toa, toacorr, Float64) - params.PLCHROMEPOCH,
        crn.ln_js,
        νref_val^params.TNCHROMIDX,
    )
end
