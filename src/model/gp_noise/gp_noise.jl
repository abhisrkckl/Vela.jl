"""
    powerlaw(A, γ, f, f1)

Powerlaw spectral model for pulsar red noise. A and γ are dimensionless,
and f and f1 are frequencies. The output has dimensions of time^2.

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

Evaluate the power law Fourier-basis red/DM/chromatic noise delay. For the latter
two, the delay corresponds to a fiducial observing frequency of 1400 MHz.
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

Evaluate the prior weights for a power law red noise. A and γ are 
dimensionless, and f1 is a frequency. Returns an array of inverse 
weights whose dimensions should be time^-2.
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

function marginalized_param_names_for_gp_noise(prefix::String, N::Int)
    # This must match the convention used in `construct_woodbury_kernel()`
    # and `pint.models.noise_model.create_fourier_design_matrix()`.
    pnames = String[]
    for fname in ("SIN", "COS")
        for ii = 1:N
            push!(pnames, "$(prefix)$(fname)_$(lpad(ii, 4, "0"))")
        end
    end
    return pnames
end
