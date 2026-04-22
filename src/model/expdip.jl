export ChromaticExponentialDip

"""
    ChromaticExponentialDip

Chromatic exponential dips (possibly due to profile change events). 
Corresponds to `SimpleExponentialDip` in `PINT`.
"""
struct ChromaticExponentialDip <: DelayComponent end

function expdip_delay(t, T, f, fref, A, γ, τ, ϵ)
    dt = t - T
    return -A * (f/fref)^γ * (τ/ϵ)^(ϵ/τ) * (τ/(τ-ϵ))^((τ-ϵ)/τ) * exp(-dt/τ) /
           (1 + exp(-dt/ϵ))
end

function delay(
    ::ChromaticExponentialDip,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)::GQ
    t = corrected_toa_value(toa, toacorr, Float64)
    f = doppler_corrected_observing_frequency(toa, toacorr)
    fref = params.EXPDIPFREF
    ϵ = params.EXPDIPEPS
    return sum(
        expdip_delay(t, T, f, fref, A, γ, τ, ϵ) for (T, A, γ, τ) in
        zip(params.EXPDIPEP_, params.EXPDIPAMP_, params.EXPDIPIDX_, params.EXPDIPTAU_)
    )
end
