function _calc_y_Ninv_y__and__logdet_N(
    inner_kernel::EcorrKernel,
    Ndiag,
    y,
    params::NamedTuple,
)
    Ntoa = length(y)
    @assert length(Ndiag) == Ntoa

    y_Ninv_y = 0.0
    logdet_Nc = 0.0
    @inbounds for group in inner_kernel.ecorr_groups
        ecorr = (group.index == 0) ? 0.0 : value(params.ECORR[group.index])
        w = ecorr * ecorr

        r_r = 0.0
        r_u = 0.0
        u_u = 0.0
        logdet_N = 0.0
        @simd for ii = group.start:group.stop
            r_r += y[ii] * y[ii] / Ndiag[ii]
            r_u += y[ii] / Ndiag[ii]
            u_u += 1 / Ndiag[ii]
            logdet_N += log(Ndiag[ii])
        end

        denom = (1 + w * u_u)
        y_Ninv_y += r_r - w * r_u * r_u / denom
        logdet_Nc += logdet_N + log(denom)
    end

    return y_Ninv_y, logdet_Nc
end

function _calc_Ninv_M(
    inner_kernel::EcorrKernel,
    M::Matrix{Float64},
    Ndiag,
    params::NamedTuple,
)
    Ntoa, Npar = size(M)
    A = Matrix{Float64}(undef, Ntoa, Npar)
    @inbounds for group in inner_kernel.ecorr_groups
        ecorr = (group.index == 0) ? 0.0 : value(params.ECORR[group.index])
        w = ecorr * ecorr
        toa_range = group.start:group.stop

        Q = 0.0
        @simd for i in toa_range
            Q += 1 / Ndiag[i]
        end

        for p = 1:Npar
            P_p = 0.0
            @simd for i in toa_range
                P_p += M[i, p] / Ndiag[i]
            end

            R = w * P_p / (1 + w * Q)

            @simd for i in toa_range
                A[i, p] = (M[i, p] - R) / Ndiag[i]
            end
        end
    end

    return A
end
