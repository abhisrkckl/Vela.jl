const fyr = frequency(1 / (365.25 * 24 * 3600))

"""The powerlaw spectrum used to model pulsar red noise.
If a_i and b_i are Fourier coefficients, the spectral
power is given by

    P(f) = <a^2> = <b^2> = A^2 / (12 * π^2 * fyr^3) * (fyr / f)^γ * f_low

The amplitude A is designed to be dimensionless for achromatic red noise.
It has dimensions of [T^-2] for DM noise. Note that this convention is 
different from the convention used in ENTERPRISE, PINT, and TEMPONEST.
There, the amplitude for DM noise is also dimensionless by construction.
"""
powerlaw(log10_A::GQ, γ::GQ, f_low::GQ, f::GQ) =
    (100.0^log10_A) / (12 * π^2 * fyr^3) * (fyr / f)^γ * f_low

ParamLocator = @NamedTuple{name::Symbol, index::Int64, default_value::Float64}
locate_param(loc::ParamLocator, params) =
    loc.index == 0 ? loc.default_value : params[loc.index]

struct PowerlawGPPrior <: Prior
    frequency::GQ{Float64}
    loc_log_amp::ParamLocator
    loc_gamma::ParamLocator
    loc_f_low::ParamLocator
end

prior_transform_param(loc, priors) =
    loc.index == 0 ? loc.default_value : prior_transform(priors[loc.index], cube[loc.index])
function prior_transform(plgp::PowerlawGPPrior, priors, cube, index::UInt)
    log10_A = prior_transform_param(plgp.loc_log_amp, priors)
    gamma = prior_transform_param(plgp.loc_gamma, priors)
    f_low = prior_transform_param(plgp.loc_f_low, priors)

    σ = powerlaw(log10_A, gamma, f_low, plgp.frequency)
    dist = Normal(0.0, σ)
    return quantile(dist, cube[index])
end

distr(plgp::PowerlawGPPrior, params) = Normal(
    0.0,
    powerlaw(
        locate_param(plgp.loc_log_amp, params),
        locate_param(plgp.loc_gamma, params),
        locate_param(plgp.loc_f_low, params),
        plgp.frequency,
    ),
)
