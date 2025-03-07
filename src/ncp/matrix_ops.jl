using LinearAlgebra
using .Threads

"""
Given the design matrix `M`, diagonal measurement noise matrix `N`, diagonal weight matrix `Φ`, and
residual vector `y`, compute the matrices
    ``Σ^-1 = Φ^-1 + M^T N^-1 M``,
and
    `u = M^T N^-1 y `.
"""
function calc_Sigmainv_and_MT_Ninv_y(
    M::Matrix{X},
    Ndiag::Vector{X},
    Phidiag::Vector{X},
    y::Vector{X},
) where {X<:AbstractFloat}
    Ntoa, Npar = size(M)
    @assert length(Ndiag) == length(y) == Ntoa
    @assert length(Phidiag) == Npar

    Ninv_M = Matrix{X}(undef, Ntoa, Npar)
    @inbounds for p = 1:Npar
        @simd for j = 1:Ntoa
            Ninv_M[j, p] = M[j, p] / Ndiag[j]
        end
    end

    # TODO: Only allocate memory for lower triangular elements.
    Sigmainv = Matrix{X}(undef, Npar, Npar)
    @inbounds for p = 1:Npar
        for q = 1:p
            Sigmainv_pq = (p == q) ? 1 / Phidiag[p] : zero(X)
            @simd for j = 1:Ntoa
                Sigmainv_pq += M[j, p] * Ninv_M[j, q]
            end

            # Only lower triangular elements are populated.
            # The rest contain garbage.
            Sigmainv[p, q] = Sigmainv_pq
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

    return Symmetric(Sigmainv, :L), u
end

"""
Given the design matrix `M`, diagonal measurement noise matrix `N`, diagonal weight matrix `Φ`, and
residual vector `y`, compute `Δa = a-a_0` such that
    ``Δa = U^-1 α -  Σ M^T N^-1 y = U^-1 (α - L^-1 M^T N^-1 y)``
where 
    ``Σ^-1 = Φ^-1 + M^T N^-1 M = L U``
and ``L = U^T``.

This function scales as O(Ntoa * Npar^2) when Ntoa >> Npar.
"""
function calc_noncentral_transform(
    M::Matrix{X},
    Ndiag::Vector{X},
    Phidiag::Vector{X},
    y::Vector{X},
    alpha::Vector{X},
) where {X<:AbstractFloat}
    Npar = length(alpha)
    @assert length(Phidiag) == Npar

    Sigmainv, MT_Ninv_y = calc_Sigmainv_and_MT_Ninv_y(M, Ndiag, Phidiag, y)

    Sigmainv_cf = cholesky!(Sigmainv)

    dα = ldiv!(Sigmainv_cf.L, MT_Ninv_y) # Linv_MT_Ninv_y
    @inbounds @simd for p = 1:Npar
        dα[p] = alpha[p] - dα[p]
    end

    return ldiv!(Sigmainv_cf.U, dα)
end
