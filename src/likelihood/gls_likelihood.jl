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

# function _calc_Σinv__MT_Ninv_y(
#     M::Matrix{X},
#     Ndiag::Vector{X},
#     Φinv::Vector{X},
#     y::Vector{X},
# ) where {X<:AbstractFloat}
#     Ntoa, Npar = size(M)
#     @assert length(Ndiag) == length(y) == Ntoa
#     @assert length(Φinv) == Npar

#     Ninv_M = Matrix{X}(undef, Ntoa, Npar)
#     @inbounds for p = 1:Npar
#         @simd for j = 1:Ntoa
#             Ninv_M[j, p] = M[j, p] / Ndiag[j]
#         end
#     end

#     # TODO: Only allocate memory for lower triangular elements.
#     Σinv = Matrix{X}(undef, Npar, Npar)
#     @inbounds for p = 1:Npar
#         for q = 1:p
#             Σinv_qp = (p == q) ? Φinv[p] : zero(X)
#             @simd for j = 1:Ntoa
#                 Σinv_qp += M[j, p] * Ninv_M[j, q]
#             end

#             # Only upper triangular elements are populated.
#             # The rest contain garbage.
#             Σinv[q, p] = Σinv_qp
#         end
#     end

#     u = Vector{X}(undef, Npar)
#     @inbounds for p = 1:Npar
#         up = zero(X)
#         @simd for j = 1:Ntoa
#             up += Ninv_M[j, p] * y[j]
#         end
#         u[p] = up
#     end

#     return Symmetric(Σinv, :U), u
# end

# function _calc_y_Ninv_y(Ndiag, y)
#     Ntoa = length(y)
#     @assert length(Ndiag) == Ntoa

#     y_Ninv_y = 0.0
#     @inbounds @simd for j in 1:Ntoa
#         y_Ninv_y += y[j] * y[j] * Ndiag[j]
#     end

#     return y_Ninv_y
# end

# function calc_lnlike_serial(
#     model::TimingModel{ComponentsTuple,WoodburyKernel{WhiteNoiseKernel,GPComponentsTuple},PriorsTuple},
#     toas::Vector{TOA},
#     params::NamedTuple
# ) where {ComponentsTuple<:Tuple, GPComponentsTuple<:Tuple, PriorsTuple<:Tuple}
#     y, Ndiag = _calc_resids_and_Ndiag(model, toas, params)
#     M = model.kernel.noise_basis
#     Φinv = calc_noise_weights_inv(model.kernel, params)

#     Σinv, MT_Ninv_y = _calc_Sigmainv__MT_Ninv_y(M, Ndiag, Phiinv, y)
#     y_Ninv_y = _calc_y_Ninv_y(Ndiag, y)

#     Σinv_cf = cholesky(Σinv)
#     Σ_MT_Ninv_y = Σinv_cf \ MT_Ninv_y
#     y_Ninv_M_Σ_MT_Ninv_y = dot(Σ_MT_Ninv_y, MT_Ninv_y)

#     logdet_Φ = -sum(log, Φinv)
#     logdet_N = sum(log, Ndiag)
#     logdet_Σinv = 2 * sum(log, diag(Σinv_cf.U))

#     return -0.5*(y_Ninv_y - y_Ninv_M_Σ_MT_Ninv_y + logdet_N + logdet_Φ + logdet_Σinv)
# end


# function calc_lnlike_serial(
#     model::TimingModel{ComponentsTuple,KernelType,PriorsTuple},
#     toas::Vector{T},
#     params::NamedTuple,
# ) where {ComponentsTuple<:Tuple,KernelType<:WoodburyKernel,PriorsTuple<:Tuple,T<:TOABase}
#     y, Ndiag = _calc_resids_and_Ndiag(model, toas, params)
#     M = model.kernel.noise_basis
#     Φ = calc_noise_weights(model, params)

#     Ninv_M = M ./ Ndiag
#     MT_Ninv_M = transpose(M) * Ninv_M
#     MT_Ninv_y = transpose(Ninv_M) * y

#     Σinv = Φ + MT_Ninv_M
#     Σinv_cf = cholesky(Σinv)

#     Σ_MT_Ninv_y = Σinv_cf \ MT_Ninv_y
#     y_Ninv_M_Σ_MT_Ninv_y = dot(MT_Ninv_y, Σ_MT_Ninv_y)

#     logdet_N = sum(log, N)
#     logdet_Φ = sum(log, Φ)
#     logdet_Σinv = 2 * sum(log, diag(a))

#     y_Ninv_y = dot(y, y ./ Ninv)

#     return -0.5 * (y_Ninv_y - y_Ninv_M_Σ_MT_Ninv_y + logdet_N + logdet_Φ + logdet_Σinv)
# end
