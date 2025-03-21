function calc_noise_weights_inv(kernel::WoodburyKernel, params::NamedTuple)
    return vcat(
        (calc_noise_weights_inv(gp_comp, params) for gp_comp in kernel.gp_components)...,
    )
end

function _calc_resids_and_Ndiag(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    ntoas = length(toas)
    ys = Vector{Float64}(undef, ntoas)
    Ndiag = Vector{Float64}(undef, ntoas)

    @inbounds for (j, toa) in enumerate(toas)
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
        ys[j] = dphase / doppler_shifted_spin_frequency(ctoa)
        Ndiag[j] = scaled_toa_error_sqr(toa, ctoa)
    end

    return ys, Ndiag
end

function _calc_resids_and_Ndiag(
    model::TimingModel,
    wtoas::Vector{WidebandTOA},
    params::NamedTuple,
)
    tzrphase = calc_tzr_phase(model, params)

    ntoas = length(wtoas)
    ys = Vector{Float64}(undef, 2 * ntoas)
    Ndiag = Vector{Float64}(undef, 2 * ntoas)

    @inbounds for (j, wtoa) in enumerate(wtoas)
        cwtoa = correct_toa(model, wtoa, params)
        dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)
        ys[j] = dphase / doppler_shifted_spin_frequency(ctoa)
        ys[ntoas+j] = dm_residual(wtoa.dminfo, cwtoa.dm_correction)
        Ndiag[j] = scaled_toa_error_sqr(wtoa.toa, cwtoa.toa_correction)
        Ndiag[ntoas+j] = scaled_dm_error_sqr(wtoa.dminfo, cwtoa.dm_correction)
    end

    return ys, Ndiag
end

function _calc_y_Ninv_y__and__logdet_N(
    ::WhiteNoiseKernel,
    Ndiag::Vector{Float64},
    y::Vector{Float64},
    ::NamedTuple,
)
    Ntoa = length(y)
    @assert length(Ndiag) == Ntoa

    y_Ninv_y = 0.0
    @inbounds @simd for j = 1:Ntoa
        y_Ninv_y += y[j] * y[j] / Ndiag[j]
    end

    logdet_N = sum(log, Ndiag)

    return y_Ninv_y, logdet_N
end

function _calc_Σinv__and__MT_Ninv_y(
    inner_kernel::Kernel,
    M::Matrix{X},
    Ndiag::Vector{X},
    Φinv::Vector{X},
    y::Vector{X},
    params::NamedTuple,
) where {X<:AbstractFloat}
    Ntoa, Npar = size(M)
    @assert length(Ndiag) == length(y) == Ntoa
    @assert length(Φinv) == Npar

    Ninv_M = _calc_Ninv_M(inner_kernel, M, Ndiag, params)

    # TODO: Only allocate memory for lower triangular elements.
    Σinv = Matrix{X}(undef, Npar, Npar)
    @inbounds for p = 1:Npar
        for q = 1:p
            Σinv_qp = (p == q) ? Φinv[p] : zero(X)
            @simd for j = 1:Ntoa
                Σinv_qp += M[j, p] * Ninv_M[j, q]
            end

            # Only upper triangular elements are populated.
            # The rest contain garbage.
            Σinv[q, p] = Σinv_qp
        end
    end

    u = Vector{X}(undef, Npar)
    @inbounds for p = 1:Npar
        up = zero(X)
        @simd for j = 1:Ntoa
            up += Ninv_M[j, p] * y[j]
        end
        u[p] = up
    end

    return Symmetric(Σinv, :U), u
end

function _calc_Ninv_M(
    ::WhiteNoiseKernel,
    M::Matrix{X},
    Ndiag::Vector{X},
    ::NamedTuple,
) where {X<:AbstractFloat}
    Ntoa, Npar = size(M)
    @assert length(Ndiag) == Ntoa

    Ninv_M = Matrix{X}(undef, Ntoa, Npar)
    @inbounds for p = 1:Npar
        @simd for j = 1:Ntoa
            Ninv_M[j, p] = M[j, p] / Ndiag[j]
        end
    end

    return Ninv_M
end

function _gls_lnlike_serial(
    inner_kernel::Kernel,
    M::Matrix{X},
    Ndiag::Vector{X},
    Φinv::Vector{X},
    y::Vector{X},
    params::NamedTuple,
) where {X<:AbstractFloat}
    Σinv, MT_Ninv_y = _calc_Σinv__and__MT_Ninv_y(inner_kernel, M, Ndiag, Φinv, y, params)
    y_Ninv_y, logdet_N = _calc_y_Ninv_y__and__logdet_N(inner_kernel, Ndiag, y, params)

    Σinv_cf = cholesky!(Σinv)
    logdet_Σinv = logdet(Σinv_cf)

    Linv_MT_Ninv_y = ldiv!(Σinv_cf.L, MT_Ninv_y)
    y_Ninv_M_Σ_MT_Ninv_y = dot(Linv_MT_Ninv_y, Linv_MT_Ninv_y)

    logdet_Φ = -sum(log, Φinv)

    return -0.5 * (y_Ninv_y - y_Ninv_M_Σ_MT_Ninv_y + logdet_N + logdet_Φ + logdet_Σinv)
end

function calc_lnlike_serial(
    model::TimingModel{
        ComponentsTuple,
        WoodburyKernel{WhiteNoiseKernel,GPComponentsTuple},
        PriorsTuple,
    },
    toas::Vector{TOAType},
    params::NamedTuple,
) where {
    ComponentsTuple<:Tuple,
    GPComponentsTuple<:Tuple,
    PriorsTuple<:Tuple,
    TOAType<:TOABase,
}
    y, Ndiag = _calc_resids_and_Ndiag(model, toas, params)
    M = model.kernel.noise_basis
    Φinv = calc_noise_weights_inv(model.kernel, params)

    return _gls_lnlike_serial(model.kernel.inner_kernel, M, Ndiag, Φinv, y, params)
end

calc_lnlike(
    model::TimingModel{
        ComponentsTuple,
        WoodburyKernel{WhiteNoiseKernel,GPComponentsTuple},
        PriorsTuple,
    },
    toas::Vector{TOAType},
    params::NamedTuple,
) where {
    ComponentsTuple<:Tuple,
    GPComponentsTuple<:Tuple,
    PriorsTuple<:Tuple,
    TOAType<:TOABase,
} = calc_lnlike_serial(model, toas, params)
