function calc_noise_weights_inv(
    kernel::WoodburyKernel,
    params_list::AbstractVector{<:NamedTuple},
)
    return map(params -> calc_noise_weights_inv(kernel, params), params_list)
end

function _calc_resids_and_Ninvdiag(
    model::TimingModel,
    toas::Vector{TOA},
    params_list::AbstractVector{<:NamedTuple},
)
    nsets = length(params_list)
    ntoas = length(toas)

    result_data = Array{Float64}(undef, ntoas, 2, nsets)

    Threads.@threads for i in eachindex(params_list)
        params = params_list[i]
        tzrphase = calc_tzr_phase(model, params)

        ys_view = @view result_data[:, 1, i]
        Ninvdiag_view = @view result_data[:, 2, i]

        @inbounds for j in eachindex(toas)
            toa = toas[j]
            ctoa = correct_toa(model, toa, params)
            dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)

            ys_view[j] = value(dphase / doppler_shifted_spin_frequency(ctoa))
            Ninvdiag_view[j] = 1.0 / value(scaled_toa_error_sqr(toa, ctoa))
        end
    end

    ys = @view result_data[:, 1, :]
    Ninvdiag = @view result_data[:, 2, :]

    return ys, Ninvdiag
end

function _calc_resids_and_Ninvdiag(
    model::TimingModel,
    wtoas::Vector{WidebandTOA},
    params_list::AbstractVector{<:NamedTuple},
)
    nsets = length(params_list)
    ntoas = length(wtoas)

    result_data = Array{Float64}(undef, 2*ntoas, 2, nsets)

    Threads.@threads for i in eachindex(params_list)
        params = params_list[i]
        tzrphase = calc_tzr_phase(model, params)

        ys_view = @view result_data[:, 1, i]
        Ninvdiag_view = @view result_data[:, 2, i]

        @inbounds for (j, wtoa) in enumerate(wtoas)
            cwtoa = correct_toa(model, wtoa, params)
            dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)

            ys_view[j] =
                value(dphase / doppler_shifted_spin_frequency(cwtoa.toa_correction))
            ys_view[ntoas+j] = value(dm_residual(wtoa.dminfo, cwtoa.dm_correction))

            Ninvdiag_view[j] =
                1.0 / value(scaled_toa_error_sqr(wtoa.toa, cwtoa.toa_correction))
            Ninvdiag_view[ntoas+j] =
                1.0 / value(scaled_dm_error_sqr(wtoa.dminfo, cwtoa.dm_correction))
        end
    end

    ys = @view result_data[:, 1, :]
    Ninvdiag = @view result_data[:, 2, :]

    return ys, Ninvdiag
end

function _calc_y_Ninv_y__and__logdet_N(
    ::WhiteNoiseKernel,
    Ninvdiag::AbstractArray,
    y::AbstractArray,
    ::AbstractVector{<:NamedTuple},
)
    @assert size(Ninvdiag) == size(y)
    Ndata, Nsets = size(y)

    result_data = Array{Float64}(undef, 2, Nsets)

    y_Ninv_y = @view result_data[1, :]
    logdet_N = @view result_data[2, :]
    for i = 1:Nsets
        y_Ninv_y_i = 0.0
        logdet_N_i = 0.0
        @inbounds @simd for j = 1:Ndata
            y_Ninv_y_i += y[j, i] * y[j, i] * Ninvdiag[j, i]
            logdet_N_i -= log(Ninvdiag[j, i])
        end # COV_EXCL_LINE

        y_Ninv_y[i] = y_Ninv_y_i
        logdet_N[i] = logdet_N_i
    end

    return y_Ninv_y, logdet_N
end

function _calc_Ninv_M(
    ::WhiteNoiseKernel,
    M::AbstractMatrix,
    Ninvdiag::AbstractMatrix,
    ::AbstractVector{<:NamedTuple};
)
    Ntoa, Npar = size(M)
    Nsets = size(Ninvdiag)[2]
    @assert size(Ninvdiag)[1] == Ntoa

    X = eltype(M)

    Ninv_M = Array{X}(undef, Ntoa, Npar, Nsets)
    Threads.@threads for i = 1:Nsets
        for p = 1:Npar
            @inbounds @simd for j = 1:Ntoa
                Ninv_M[j, p, i] = M[j, p] * Ninvdiag[j, i]
            end # COV_EXCL_LINE
        end
    end

    return Ninv_M
end












function _calc_Σinv__and__MT_Ninv_y(
    inner_kernel::Kernel,
    M::AbstractMatrix,
    Ninvdiag::AbstractMatrix,
    Φinv::AbstractVector,
    y::AbstractVector,
    params::NamedTuple,
)
    Ntoa, Npar = size(M)
    @assert length(Ninvdiag) == length(y) == Ntoa
    @assert length(Φinv) == Npar

    Ninv_M = _calc_Ninv_M(inner_kernel, M, Ninvdiag, params)

    X = eltype(M)

    # TODO: Only allocate memory for lower triangular elements.
    result_data = Matrix{X}(undef, Npar, Npar + 1)
    Σinv = @view result_data[:, 1:Npar]
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

    u = @view result_data[:, Npar+1]
    @inbounds @threads for p = 1:Npar
        up = zero(X)
        @simd for j = 1:Ntoa
            up += Ninv_M[j, p] * y[j]
        end # COV_EXCL_LINE
        u[p] = up
    end # COV_EXCL_LINE

    return Symmetric(Σinv, :U), u
end
