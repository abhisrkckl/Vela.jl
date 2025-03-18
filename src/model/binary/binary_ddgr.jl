export BinaryDDGR

"""The incarnation of the Damour & Deruelle binary model for eccentric binaries assuming
general relativity.

Reference:
    [Taylor & Weisberg 1989](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..434T/abstract)
"""
struct BinaryDDGR <: BinaryDDBase end

"""Setting this function to zero corresponds to the 1PN-accurate version of
Kepler's third law. Also returns the first derivative."""
function _ar_func_and_deriv(ar, M, m2, n)
    m1 = M - m2
    η = m1 * m2 / (M * M)

    f_lhs = ar * ar * ar
    f_rhs_0PN = M / (n * n)
    f_rhs_1PN_fac = 1 + (η - 9) * M / (2 * ar)
    f = f_lhs - f_rhs_0PN * f_rhs_1PN_fac * f_rhs_1PN_fac

    fp_lhs = 3 * ar * ar
    fp_rhs_fac = (η - 9) * M / (-2 * ar * ar)
    fp = fp_lhs - f_rhs_0PN * 2 * f_rhs_1PN_fac * fp_rhs_fac

    return f, fp
end

"""Solve the 1PN-accurate Kepler's 3rd law using Newton-Raphson method."""
function _kepler_3rd_law_1pn(M, m2, n, ar_init)
    ar0 = distance(-1.0)
    ar1 = ar_init

    while abs((ar1 - ar0) / ar1) > 1e-10
        f, fp = _ar_func_and_deriv(ar1, M, m2, n)
        ar0 = ar1
        ar1 = ar0 - (f / fp)
    end

    return ar1
end

function mean_motion_ddgr()
    
end

function DDState(ddgr::BinaryDDGR, toa::TOA, toacorr::TOACorrection, params::NamedTuple)
    Δt = corrected_toa_value(toa, toacorr, Float64) - params.T0
    n = mean_motion(Δt, params, dd.use_fbx)
    l = mean_anomaly(Δt, params, dd.use_fbx)

    M = params.MTOT
    m2 = params.m2

    a1 = params.A1 + Δt * params.A1DOT
    et = params.ECC + Δt * params.EDOT

    ar0 = cbrt(M / (n * n))
    ar = _kepler_3rd_law_1pn(M, m2, n, ar0)
    ap = ar * m2 / M
    sini = a1 / ap

    # This expression uses ar0 instead of ar to be consistent with tempo and PINT.
    # tempo2 uses ar in this context. In the tempo source code this difference is 
    # attributed to "advice from Norbert [Wex]".
    γ = et * m2 * (M + m2) / (n * ar0 * M)


end