export calc_lnlike_vectorized

function calc_noise_weights_inv(kernel::WoodburyKernel, params::NamedTuple)
    return vcat(
        (calc_noise_weights_inv(gp_comp, params) for gp_comp in kernel.gp_components)...,
    )
end

function _calc_resids_and_Ninvdiag(
    model::TimingModel,
    toas::Vector{TOA},
    params::NamedTuple;
    parallel::Bool = false,
)
    tzrphase = calc_tzr_phase(model, params)

    ntoas = length(toas)

    result_data = Vector{Float64}(undef, 2 * ntoas)
    ys = @view result_data[1:ntoas]
    Ninvdiag = @view result_data[(ntoas+1):end]

    if parallel
        @inbounds @threads for j = 1:ntoas
            toa = toas[j]
            ctoa = correct_toa(model, toa, params)
            dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
            ys[j] = value(dphase / doppler_shifted_spin_frequency(ctoa))
            Ninvdiag[j] = 1.0 / value(scaled_toa_error_sqr(toa, ctoa))
        end # COV_EXCL_LINE
    else
        @inbounds for j = 1:ntoas
            toa = toas[j]
            ctoa = correct_toa(model, toa, params)
            dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
            ys[j] = value(dphase / doppler_shifted_spin_frequency(ctoa))
            Ninvdiag[j] = 1.0 / value(scaled_toa_error_sqr(toa, ctoa))
        end
    end

    return ys, Ninvdiag
end

function _calc_resids_and_Ninvdiag(
    model::TimingModel,
    wtoas::Vector{WidebandTOA},
    params::NamedTuple;
    parallel::Bool = false,
)
    tzrphase = calc_tzr_phase(model, params)

    ntoas = length(wtoas)
    result_data = Vector{Float64}(undef, 4 * ntoas)
    ys = @view result_data[1:(2*ntoas)]
    Ninvdiag = @view result_data[(2*ntoas+1):end]

    if parallel
        @inbounds @threads for j = 1:ntoas
            wtoa = wtoas[j]
            cwtoa = correct_toa(model, wtoa, params)
            dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)
            ys[j] = dphase / doppler_shifted_spin_frequency(cwtoa.toa_correction)
            ys[ntoas+j] = dm_residual(wtoa.dminfo, cwtoa.dm_correction)
            Ninvdiag[j] = 1.0 / scaled_toa_error_sqr(wtoa.toa, cwtoa.toa_correction)
            Ninvdiag[ntoas+j] = 1.0 / scaled_dm_error_sqr(wtoa.dminfo, cwtoa.dm_correction)
        end # COV_EXCL_LINE
    else
        @inbounds for (j, wtoa) in enumerate(wtoas)
            cwtoa = correct_toa(model, wtoa, params)
            dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)
            ys[j] = dphase / doppler_shifted_spin_frequency(cwtoa.toa_correction)
            ys[ntoas+j] = dm_residual(wtoa.dminfo, cwtoa.dm_correction)
            Ninvdiag[j] = 1.0 / scaled_toa_error_sqr(wtoa.toa, cwtoa.toa_correction)
            Ninvdiag[ntoas+j] = 1.0 / scaled_dm_error_sqr(wtoa.dminfo, cwtoa.dm_correction)
        end
    end

    return ys, Ninvdiag
end

function _calc_y_Ninv_y__and__logdet_N(
    ::WhiteNoiseKernel,
    Ninvdiag::AbstractVector,
    y::AbstractVector,
    ::NamedTuple;
    parallel::Bool = false,
)
    Ntoa = length(y)
    @assert length(Ninvdiag) == Ntoa

    y_Ninv_y = 0.0
    @inbounds @simd for j = 1:Ntoa
        y_Ninv_y += y[j] * y[j] * Ninvdiag[j]
    end # COV_EXCL_LINE

    logdet_N = -sum(log, Ninvdiag)

    return y_Ninv_y, logdet_N
end

function _calc_Σinv__and__MT_Ninv_y(
    inner_kernel::Kernel,
    M::AbstractMatrix,
    Ninvdiag::AbstractVector,
    Φinv::AbstractVector,
    y::AbstractVector,
    params::NamedTuple;
    parallel::Bool = false,
)
    Ntoa, Npar = size(M)
    @assert length(Ninvdiag) == length(y) == Ntoa
    @assert length(Φinv) == Npar

    Ninv_M = _calc_Ninv_M(inner_kernel, M, Ninvdiag, params)

    X = eltype(M)

    # TODO: Only allocate memory for lower triangular elements.
    result_data = Matrix{X}(undef, Npar, Npar + 1)
    Σinv = @view result_data[:, 1:Npar]
    if parallel
        @inbounds @threads for p = 1:Npar
            for q = 1:p
                Σinv_qp = (p == q) ? Φinv[p] : zero(Φinv[p])
                @simd for j = 1:Ntoa
                    Σinv_qp += M[j, p] * Ninv_M[j, q]
                end # COV_EXCL_LINE

                # Only upper triangular elements are populated.
                # The rest contain garbage.
                Σinv[q, p] = Σinv_qp
            end
        end # COV_EXCL_LINE
    else
        @inbounds for p = 1:Npar
            for q = 1:p
                Σinv_qp = (p == q) ? Φinv[p] : zero(Φinv[p])
                @simd for j = 1:Ntoa
                    Σinv_qp += M[j, p] * Ninv_M[j, q]
                end # COV_EXCL_LINE

                # Only upper triangular elements are populated.
                # The rest contain garbage.
                Σinv[q, p] = Σinv_qp
            end
        end
    end

    u = @view result_data[:, Npar+1]
    if parallel
        @inbounds @threads for p = 1:Npar
            up = zero(X)
            @simd for j = 1:Ntoa
                up += Ninv_M[j, p] * y[j]
            end # COV_EXCL_LINE
            u[p] = up
        end # COV_EXCL_LINE
    else
        @inbounds for p = 1:Npar
            up = zero(X)
            @simd for j = 1:Ntoa
                up += Ninv_M[j, p] * y[j]
            end # COV_EXCL_LINE
            u[p] = up
        end
    end

    return Symmetric(Σinv, :U), u
end

function _calc_Ninv_M(
    ::WhiteNoiseKernel,
    M::AbstractMatrix,
    Ninvdiag::AbstractVector,
    ::NamedTuple;
    parallel::Bool = false,
)
    Ntoa, Npar = size(M)
    @assert length(Ninvdiag) == Ntoa

    X = eltype(M)

    Ninv_M = Matrix{X}(undef, Ntoa, Npar)
    if parallel
        @inbounds @threads for p = 1:Npar
            @simd for j = 1:Ntoa
                Ninv_M[j, p] = M[j, p] * Ninvdiag[j]
            end # COV_EXCL_LINEi
        end # COV_EXCL_LINE
    else
        @inbounds for p = 1:Npar
            @simd for j = 1:Ntoa
                Ninv_M[j, p] = M[j, p] * Ninvdiag[j]
            end # COV_EXCL_LINE
        end
    end

    return Ninv_M
end

function _gls_lnlike_serial(
    inner_kernel::Kernel,
    M::AbstractMatrix,
    Ninvdiag::AbstractVector,
    Φinv::AbstractVector,
    y::AbstractVector,
    params::NamedTuple,
)
    Σinv, MT_Ninv_y = _calc_Σinv__and__MT_Ninv_y(inner_kernel, M, Ninvdiag, Φinv, y, params)

    if !isposdef(Σinv)
        return -Inf # COV_EXCL_LINE
    end

    y_Ninv_y, logdet_N = _calc_y_Ninv_y__and__logdet_N(inner_kernel, Ninvdiag, y, params)

    Σinv_cf = cholesky!(Σinv)
    logdet_Σinv = logdet(Σinv_cf)

    Linv_MT_Ninv_y = ldiv!(Σinv_cf.L, MT_Ninv_y)
    y_Ninv_M_Σ_MT_Ninv_y = dot(Linv_MT_Ninv_y, Linv_MT_Ninv_y)

    logdet_Φ = -sum(log, Φinv)

    return -0.5 * (y_Ninv_y - y_Ninv_M_Σ_MT_Ninv_y + logdet_N + logdet_Φ + logdet_Σinv)
end

function _gls_lnlike_parallel(
    inner_kernel::Kernel,
    M::AbstractMatrix,
    Ninvdiag::AbstractVector,
    Φinv::AbstractVector,
    y::AbstractVector,
    params::NamedTuple,
)
    Σinv, MT_Ninv_y = _calc_Σinv__and__MT_Ninv_y(
        inner_kernel,
        M,
        Ninvdiag,
        Φinv,
        y,
        params;
        parallel = true,
    )

    if !isposdef(Σinv)
        return -Inf # COV_EXCL_LINE
    end

    y_Ninv_y, logdet_N =
        _calc_y_Ninv_y__and__logdet_N(inner_kernel, Ninvdiag, y, params; parallel = true)

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
    y, Ninvdiag = _calc_resids_and_Ninvdiag(model, toas, params)
    M = model.kernel.noise_basis
    Φinv = calc_noise_weights_inv(model.kernel, params)

    return _gls_lnlike_serial(model.kernel.inner_kernel, M, Ninvdiag, Φinv, y, params)
end

function calc_lnlike(
    model::TimingModel{ComponentsTuple,WoodburyKernelType,PriorsTuple},
    toas::Vector{TOAType},
    params::NamedTuple,
) where {
    ComponentsTuple<:Tuple,
    WoodburyKernelType<:WoodburyKernel,
    PriorsTuple<:Tuple,
    TOAType<:TOABase,
}
    y, Ninvdiag = _calc_resids_and_Ninvdiag(model, toas, params; parallel = true)
    M = model.kernel.noise_basis
    Φinv = calc_noise_weights_inv(model.kernel, params)

    return _gls_lnlike_parallel(model.kernel.inner_kernel, M, Ninvdiag, Φinv, y, params)
end

function calc_lnlike_vectorized(
    model::TimingModel,
    toas::Vector{T},
    paramss,
) where {T<:TOABase}
    nparamss = size(paramss)[1]
    result = Vector{Float64}(undef, nparamss)
    @threads :static for ii = 1:nparamss
        result[ii] = calc_lnlike_serial(model, toas, paramss[ii, :])
    end # COV_EXCL_LINE
    return result
end
