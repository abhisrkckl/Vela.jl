function calc_noise_weights_inv(kernel::WoodburyKernel, params::NamedTuple)
    return vcat(
        (calc_noise_weights_inv(gp_comp, params) for gp_comp in kernel.gp_components)...,
    )
end


# function calc_noise_weights_inv(kernel::Vela.WoodburyKernel, params::NamedTuple)
#     Np = size(kernel.noise_basis)[2]
#     weights_inv = Vector{Float64}(undef, Np)

#     start_idx = 1
#     for gp_comp in kernel.gp_components
#         start_idx += Vela.calc_noise_weights_inv!(gp_comp, params, start_idx, weights_inv)
#     end

#     return weights_inv
# end

# function _calc_resids_and_Ndiag(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
#     tzrphase = calc_tzr_phase(model, params)

#     ntoas = length(toas)
#     ys = Vector{Float64}(undef, ntoas)
#     Ndiag = Vector{Float64}(undef, ntoas)

#     @inbounds for (j, toa) in enumerate(toas)
#         ctoa = correct_toa(model, toa, params)
#         dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
#         ys[j] = dphase / doppler_shifted_spin_frequency(ctoa)
#         Ndiag[j] = scaled_toa_error_sqr(toa, ctoa)
#     end

#     return ys, Ndiag
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
