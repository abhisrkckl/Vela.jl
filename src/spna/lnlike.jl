function calc_y_and_Ninvdiag(spna::NarrowbandSPNA, params::NamedTuple)
    n = length(spna.resids)

    result_data = Vector{Float64}(undef, 2 * n)
    ys = @view result_data[1:n]
    Ninvdiag = @view result_data[(n+1):end]

    for ii = 1:n
        res = spna.resids[ii]
        cres = correct_resid(spna.components, res, NarrowbandResidCorrection(), params, ii)
        ys[ii] = corrected_time_residual(res, cres)
        Ninvdiag[ii] = 1 / scaled_toa_error_sqr(res, cres)
    end

    return ys, Ninvdiag
end

function calc_y_and_Ninvdiag(spna::WidebandSPNA, params::NamedTuple)
    n = length(spna.resids)

    result_data = Vector{Float64}(undef, 4 * n)
    rs = @view result_data[1:n]
    ds = @view result_data[(n+1):(2*n)]
    ςs = @view result_data[(2*n+1):(3*n)]
    ϵs = @view result_data[(3*n+1):(4*n)]

    ys = @view result_data[1:(2*n)]
    Ninvdiag = @view result_data[(2*n+1):end]

    for ii = 1:n
        res = spna.resids[ii]
        cres = correct_resid(spna.components, res, WidebandResidCorrection(), params, ii)
        rs[ii] = corrected_time_residual(res, cres)
        ds[ii] = corrected_dm_residual(res, cres)
        ςs[ii] = 1 / scaled_toa_error_sqr(res, cres)
        ϵs[ii] = 1 / scaled_dm_error_sqr(res, cres)
    end

    return ys, Ninvdiag
end

function calc_lnlike_serial(spna::SPNABase, params::NamedTuple)
    y, Ninvdiag = calc_y_and_Ninvdiag(spna, params)
    M = spna.kernel.noise_basis
    Φinv = calc_noise_weights_inv(spna.kernel, params)

    return _gls_lnlike_serial(model.kernel.inner_kernel, M, Ninvdiag, Φinv, y, params)
end
