function calc_Sigmainv_and_MT_Cinv_y__ecorr(
    M::Matrix{X},
    Ndiag::Vector{X},
    Phidiag::Vector{X},
    y::Vector{X},
    cs::Vector{X},
    groups::Vector{EcorrGroup},
) where {X<:AbstractFloat}
    Ntoa, Npar = size(M)
    @assert length(Ndiag) == length(y) == Ntoa
    @assert length(Phidiag) == Npar
    @assert length(cs) == length(groups)

    # Only upper diagonal elements are set. Rest are garbage.
    Sigmainv = Matrix{X}(undef, Npar, Npar)
    @inbounds for p in 1:Npar
        @simd for q in 1:(p-1)
            Sigmainv[q,p] = zero(X)
        end
        Sigmainv[p,p] = 1 / Phidiag[p]
    end

    u = zeros(Npar)

    @inbounds for (cg, group) in zip(cs, groups)
        toa_range = group.start:group.stop
        
        r_g = zero(X)
        @simd for j in toa_range
            r_g += 1 / Ndiag[i]
        end
        r_g *= cg * cg

        q_g = zeros(X, )

        for p in 1:Npar
            q_p = zero(X)
            for q in 1:Npar
                Sigmainv_g_qp = zero(X)
                P_qp = zero(X)
                @simd for j in group.start:group.stop

                end
                Sigmainv[q,p] += Sigmainv_g_qp
            end
        end
    end
end