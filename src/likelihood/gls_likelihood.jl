function calc_noise_weights_inv(kernel::WoodburyKernel, params::NamedTuple)
    return vcat(
        (calc_noise_weights_inv(gp_comp, params) for gp_comp in kernel.gp_components)...,
    )
end

function _calc_resids_and_Ndiag(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    ntoas = length(toas)

    result_data = Vector{Float64}(undef, 2 * ntoas)
    ys = @view result_data[1:ntoas]
    Ndiag = @view result_data[(ntoas+1):end]

    @inbounds for j = 1:ntoas
        toa = toas[j]
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
        ys[j] = value(dphase / doppler_shifted_spin_frequency(ctoa))
        Ndiag[j] = value(scaled_toa_error_sqr(toa, ctoa))
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
    result_data = Vector{Float64}(undef, 4 * ntoas)
    ys = @view result_data[1:(2*ntoas)]
    Ndiag = @view result_data[(2*ntoas+1):end]

    @inbounds for (j, wtoa) in enumerate(wtoas)
        cwtoa = correct_toa(model, wtoa, params)
        dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)
        ys[j] = dphase / doppler_shifted_spin_frequency(cwtoa.toa_correction)
        ys[ntoas+j] = dm_residual(wtoa.dminfo, cwtoa.dm_correction)
        Ndiag[j] = scaled_toa_error_sqr(wtoa.toa, cwtoa.toa_correction)
        Ndiag[ntoas+j] = scaled_dm_error_sqr(wtoa.dminfo, cwtoa.dm_correction)
    end

    return ys, Ndiag
end

function _calc_y_Ninv_y__and__logdet_N(
    ::WhiteNoiseKernel,
    Ndiag::AbstractVector,
    y::AbstractVector,
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
    M::AbstractMatrix,
    Ndiag::AbstractVector,
    Φinv::AbstractVector,
    y::AbstractVector,
    params::NamedTuple,
)
    Ntoa, Npar = size(M)
    @assert length(Ndiag) == length(y) == Ntoa
    @assert length(Φinv) == Npar

    Ninv_M = _calc_Ninv_M(inner_kernel, M, Ndiag, params)

    X = eltype(M)

    # TODO: Only allocate memory for lower triangular elements.
    result_data = Matrix{X}(undef, Npar, Npar + 1)
    Σinv = @view result_data[:, 1:Npar]
    @inbounds for p = 1:Npar
        for q = 1:p
            Σinv_qp = (p == q) ? Φinv[p] : zero(Φinv[p])
            @simd for j = 1:Ntoa
                Σinv_qp += M[j, p] * Ninv_M[j, q]
            end

            # Only upper triangular elements are populated.
            # The rest contain garbage.
            Σinv[q, p] = Σinv_qp
        end
    end

    u = @view result_data[:, Npar+1]
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
    M::AbstractMatrix,
    Ndiag::AbstractVector,
    ::NamedTuple,
)
    Ntoa, Npar = size(M)
    @assert length(Ndiag) == Ntoa

    X = eltype(M)

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
    M::AbstractMatrix,
    Ndiag::AbstractVector,
    Φinv::AbstractVector,
    y::AbstractVector,
    params::NamedTuple,
)
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
    model::TimingModel{ComponentsTuple,WoodburyKernelType,PriorsTuple},
    toas::Vector{TOAType},
    params::NamedTuple,
) where {
    ComponentsTuple<:Tuple,
    WoodburyKernelType<:WoodburyKernel,
    PriorsTuple<:Tuple,
    TOAType<:TOABase,
}
    y, Ndiag = _calc_resids_and_Ndiag(model, toas, params)
    M = model.kernel.noise_basis
    Φinv = calc_noise_weights_inv(model.kernel, params)

    return _gls_lnlike_serial(model.kernel.inner_kernel, M, Ndiag, Φinv, y, params)
end

calc_lnlike(
    model::TimingModel{ComponentsTuple,WoodburyKernelType,PriorsTuple},
    toas::Vector{TOAType},
    params::NamedTuple,
) where {
    ComponentsTuple<:Tuple,
    WoodburyKernelType<:WoodburyKernel,
    PriorsTuple<:Tuple,
    TOAType<:TOABase,
} = calc_lnlike_serial(model, toas, params)
