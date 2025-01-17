abstract type KernelSection end

struct MarginalizedPowerlawRedNoiseGP{N} <: KernelSection
    basis::Matrix{Float64}
    ln_js::NTuple{N,Float32}
end

function weight(rn::MarginalizedPowerlawRedNoiseGP, params)
    log10_A = params.TNREDAMP
    γ = params.TNREDGAM
    A = exp10(log10_A)
    f1 = params.PLREDFREQ
    w1 = powerlaw(A, γ, rn.f1, rn.f1)
    return map(ln_j -> w1 * exp(-γ * ln_j), rn.ln_js)
end

struct WoodburyKernel{InnerKernel<:Kernel, SectionsTuple<:Tuple} <: Kernel
    inner_kernel::InnerKernel
    sections::SectionsTuple
    basis::Matrix{Float64}

    function WoodburyKernel(inner_kernel::Kernel, sections::Tuple)
        @assert all(map(s -> isa(s, KernelSection), sections))
        basis = combine_bases(sections)
        return new(inner_kernel, sections, basis)
    end
end

function combine_bases(sections)
    return hcat(map(section -> section.basis, sections)...)
end

function combine_weights(sections, params)
    
end

function apply_inner_kernel(
    ::WhiteNoiseKernel,
    model::TimingModel,
    toas::Vector{TOA},
    params::NamedTuple,
    tzrphase::GQ,
    basis::Matrix{Float64},
    weights::Matrix{GQ{2,Float64}},
)
    ntoa = length(toas)
    nbasis = size(basis)[2]
    @assert size(basis)[1] == ntoa

    logdet_N = 0.0
    r_Ninv_r = 0.0
    UT_Ninv_r = zeros(GQ{-1,Float64}, nbasis)
    Σ = zeros(GQ{-2,Float64}, nbasis, nbasis)

    for (ii, toa) in enumerate(toas)
        ctoa = correct_toa(model, toa, params)

        ς2 = scaled_toa_error_sqr(toa, ctoa)
        logdet_N += log(value(ς2))

        dϕ = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
        r = dϕ / doppler_shifted_spin_frequency(ctoa)
        r_Ninv_r += value(r * r / ς2)

        for jj = 1:nbasis
            UT_Ninv_r[jj] += basis[ii, jj] * r / ς2

            for kk = 1:jj
                Σ[jj, kk] += basis[ii, jj] * basis[ii, kk] / ς2
            end
        end
    end

    for jj = 1:nbasis
        Σ[jj, jj] += 1 / weights[jj]

        for kk = (jj+1):nbasis
            Σ[jj, kk] += Σ[kk, jj]
        end
    end

    return logdet_N, r_Ninv_r, UT_Ninv_r, Σ
end
