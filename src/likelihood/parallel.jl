function calc_lnpost_vectorized(
    model::TimingModel,
    toas::Vector{T},
    paramss,
) where {T<:TOABase}
    nparamss = size(paramss, 1)

    lnposts = Vector{Float64}(undef, nparamss)

    # == Pre-allocating memory ==================================
    nthr = nthreads()

    ndata = get_ndata(toas)
    nmpar = get_nmpar(model.kernel)

    yN_size = 2*ndata
    Ninv_M_size = ndata*nmpar
    Sigmainv_size = nmpar*nmpar
    MT_Ninv_y_size = nmpar

    memory_buffer = Array{Float64}(
        undef,
        (yN_size + Ninv_M_size + Sigmainv_size + MT_Ninv_y_size),
        nthr,
    )

    y_Ninvdiag_idx = 1:yN_size
    Ninv_M_idx = (yN_size+1):(yN_size+Ninv_M_size)
    Sigmainv_idx = (yN_size+Ninv_M_size+1):(yN_size+Ninv_M_size+Sigmainv_size)
    MT_Ninv_y_idx =
        (yN_size+Ninv_M_size+Sigmainv_size+1):(yN_size+Ninv_M_size+Sigmainv_size+MT_Ninv_y_size)
    # ===========================================================

    M = model.kernel.noise_basis

    @threads for ii = 1:nparamss
        params = read_params(model, paramss[ii, :])

        lnpr = calc_lnprior(model, params)
        if !isfinite(lnpr)
            lnposts[ii] = -Inf
            continue
        end

        ithr = findfirst(
            x -> x==threadid(),
            Base.Threads.threadpooltids(Base.Threads.threadpool()),
        )

        y_Ninvdiag_buf = @view memory_buffer[y_Ninvdiag_idx, ithr]
        Ninv_M_buf = @view memory_buffer[Ninv_M_idx, ithr]
        Σinv_buf = @view memory_buffer[Sigmainv_idx, ithr]
        MT_Ninv_y_buf = @view memory_buffer[MT_Ninv_y_idx, ithr]

        Φinv = calc_noise_weights_inv(model.kernel, params)

        try
            y, Ninvdiag = _calc_resids_and_Ninvdiag!(model, toas, params, y_Ninvdiag_buf)

            Ninv_M =
                _calc_Ninv_M!(model.kernel.inner_kernel, M, Ninvdiag, params, Ninv_M_buf)

            Σinv = _calc_Σinv!(M, Ninv_M, Φinv, y, Σinv_buf)

            MT_Ninv_y = _calc_MT_Ninv_y!(Ninv_M, y, MT_Ninv_y_buf)

            y_Ninv_y, logdet_N = _calc_y_Ninv_y__and__logdet_N(
                model.kernel.inner_kernel,
                Ninvdiag,
                y,
                params,
            )

            Σinv_cf = cholesky!(Σinv)

            logdet_Σinv = logdet(Σinv_cf)

            Linv_MT_Ninv_y = ldiv!(Σinv_cf.L, MT_Ninv_y)
            y_Ninv_M_Σ_MT_Ninv_y = dot(Linv_MT_Ninv_y, Linv_MT_Ninv_y)

            logdet_Φ = -sum(log, Φinv)

            lnl =
                -0.5 * (y_Ninv_y - y_Ninv_M_Σ_MT_Ninv_y + logdet_N + logdet_Φ + logdet_Σinv)

            lnposts[ii] = lnpr + lnl
        catch
            lnposts[ii] = -Inf
            continue
        end
    end # COV_EXCL_LINE

    return lnposts
end

get_ndata(toas::Vector{TOA}) = length(toas)
get_ndata(toas::Vector{WidebandTOA}) = 2*length(toas)

get_nmpar(kernel::WoodburyKernel) = size(kernel.noise_basis, 2)

function _calc_resids_and_Ninvdiag!(
    model::TimingModel,
    toas::Vector{TOA},
    params::NamedTuple,
    y_Ninvdiag_buf::AbstractVector,
)
    tzrphase = calc_tzr_phase(model, params)

    ntoas = length(toas)

    ys = @view y_Ninvdiag_buf[1:ntoas]
    Ninvdiag = @view y_Ninvdiag_buf[(ntoas+1):end]

    @inbounds for j = 1:ntoas
        toa = toas[j]
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
        ys[j] = value(dphase / doppler_shifted_spin_frequency(ctoa))
        Ninvdiag[j] = 1.0 / value(scaled_toa_error_sqr(toa, ctoa))
    end

    return ys, Ninvdiag
end

function _calc_resids_and_Ninvdiag!(
    model::TimingModel,
    wtoas::Vector{WidebandTOA},
    params::NamedTuple,
    y_Ninvdiag_buf::AbstractVector,
)
    tzrphase = calc_tzr_phase(model, params)

    ntoas = length(wtoas)
    ys = @view y_Ninvdiag_buf[1:(2*ntoas)]
    Ninvdiag = @view y_Ninvdiag_buf[(2*ntoas+1):end]

    @inbounds for (j, wtoa) in enumerate(wtoas)
        cwtoa = correct_toa(model, wtoa, params)
        dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)
        ys[j] = value(dphase / doppler_shifted_spin_frequency(cwtoa.toa_correction))
        ys[ntoas+j] = value(dm_residual(wtoa.dminfo, cwtoa.dm_correction))
        Ninvdiag[j] = 1.0 / value(scaled_toa_error_sqr(wtoa.toa, cwtoa.toa_correction))
        Ninvdiag[ntoas+j] =
            1.0 / value(scaled_dm_error_sqr(wtoa.dminfo, cwtoa.dm_correction))
    end

    return ys, Ninvdiag
end

function _calc_Ninv_M!(
    ::WhiteNoiseKernel,
    M::AbstractMatrix,
    Ninvdiag::AbstractVector,
    ::NamedTuple,
    Ninv_M_buf::AbstractVector,
)
    Ntoa, Npar = size(M)
    @assert length(Ninvdiag) == Ntoa

    Ninv_M = reshape(Ninv_M_buf, Ntoa, Npar)
    @inbounds for p = 1:Npar
        @simd for j = 1:Ntoa
            Ninv_M[j, p] = M[j, p] * Ninvdiag[j]
        end # COV_EXCL_LINE
    end

    return Ninv_M
end

function _calc_Ninv_M!(
    inner_kernel::EcorrKernel,
    M::Matrix{Float64},
    Ninvdiag,
    params::NamedTuple,
    Ninv_M_buf::AbstractVector,
)
    Ntoa, Npar = size(M)
    A = reshape(Ninv_M_buf, Ntoa, Npar)

    @inbounds for group in inner_kernel.ecorr_groups
        ecorr = (group.index == 0) ? 0.0 : value(params.ECORR[group.index])
        w = ecorr * ecorr
        toa_range = group.start:group.stop

        Q = 0.0
        @simd for i in toa_range
            Q += Ninvdiag[i]
        end # COV_EXCL_LINE

        α = w / (1 + w * Q)

        for p = 1:Npar
            P_p = 0.0
            @simd for i in toa_range
                P_p += M[i, p] * Ninvdiag[i]
            end # COV_EXCL_LINE

            R = P_p * α

            @simd for i in toa_range
                A[i, p] = (M[i, p] - R) * Ninvdiag[i]
            end # COV_EXCL_LINE
        end
    end

    return A
end


function _calc_Σinv!(
    M::AbstractMatrix,
    Ninv_M::AbstractMatrix,
    Φinv::AbstractVector,
    y::AbstractVector,
    Σinv_buf::AbstractVector,
)
    Ntoa, Npar = size(Ninv_M)
    @assert length(y) == Ntoa
    @assert length(Φinv) == Npar
    @assert size(M) == (Ntoa, Npar)

    # TODO: Only allocate memory for lower triangular elements.
    Σinv = reshape(Σinv_buf, Npar, Npar)
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

    return Symmetric(Σinv, :U)
end

function _calc_MT_Ninv_y!(
    Ninv_M::AbstractMatrix,
    y::AbstractVector,
    MT_Ninv_y_buf::AbstractVector,
)
    Ntoa, Npar = size(Ninv_M)
    @assert length(y) == Ntoa

    u = reshape(MT_Ninv_y_buf, Npar)
    @inbounds for p = 1:Npar
        up = 0.0
        @simd for j = 1:Ntoa
            up += Ninv_M[j, p] * y[j]
        end # COV_EXCL_LINE
        u[p] = up
    end

    return u
end
