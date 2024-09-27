export PowerlawRedNoiseGP

struct PowerlawRedNoiseGP <: DelayComponent end

function powerlaw(A, γ, f, f1)
    fyr = frequency(1 / 3600 / 24 / 365.25)
    denom = 12 * π^2 * fyr * fyr * fyr
    return A * A / denom * (fyr/f)^γ * f1
end

function delay(::PowerlawRedNoiseGP, toa::TOA, toacorr::TOACorrection, params::NamedTuple)
    αs = params.PLREDSIN_
    βs = params.PLREDCOS_
    f1 = params.PLREDFREQ
    t0 = params.PLREDEPOCH

    @assert length(αs) == length(βs)

    A = 10^params.TNREDAMP
    γ = params.TNREDGAM

    Δt = corrected_toa_value(toa, toacorr, Float64) - t0

    result = time(0.0)
    for (ii, (α, β)) in enumerate(zip(αs, βs))
        f = ii * f1
        σ = sqrt(powerlaw(A, γ, f, f1))
        a = σ * α
        b = σ * β
        ϕ = 2π * f * Δt
        result += dot((a, b), sincos(ϕ))
    end

    return result
end

