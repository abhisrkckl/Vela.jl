function calc_noise_weights_inv(
    kernel::WoodburyKernel,
    params_list::AbstractVector{<:NamedTuple},
)
    return tuple(
        (
            calc_noise_weights_inv(gp_comp, params) for params in params_list for
            gp_comp in kernel.gp_components
        )...,
    )
end

function calc_resids_and_Ninvdiag(
    model::TimingModel,
    toas::Vector{TOA},
    params_list::AbstractVector{<:NamedTuple},
)
    nsets = length(params_list)
    ntoas = length(toas)

    result_data = Vector{Float64}(undef, 2 * nsets * ntoas)
    ys = @view result_data[1:(nsets*ntoas)]
    Ninvdiag = @view result_data[(nsets*ntoas+1):end]

    Threads.@threads for i = 1:nsets
        params = params_list[i]
        tzrphase = calc_tzr_phase(model, params)

        inds = ((i-1)*ntoas+1):(i*ntoas)
        ys_view = @view ys[inds]
        Ninvdiag_view = @view Ninvdiag[inds]

        @inbounds for (j, toa) in enumerate(toas)
            ctoa = correct_toa(model, toa, params)
            dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)

            ys_view[j] = value(dphase / doppler_shifted_spin_frequency(ctoa))
            Ninvdiag_view[j] = 1.0 / value(scaled_toa_error_sqr(toa, ctoa))
        end
    end

    return ys, Ninvdiag
end

function calc_resids_and_Ninvdiag(
    model::TimingModel,
    wtoas::Vector{WidebandTOA},
    params_list::AbstractVector{<:NamedTuple},
)
    nsets = length(params_list)
    ntoas = length(wtoas)

    result_data = Vector{Float64}(undef, 4 * nsets * ntoas)
    ys = @view result_data[1:(2*nsets*ntoas)]
    Ninvdiag = @view result_data[(2*nsets*ntoas+1):end]

    Threads.@threads for i = 1:nsets
        params = params_list[i]
        tzrphase = calc_tzr_phase(model, params)

        inds = ((i-1)*2*ntoas+1):(i*2*ntoas)
        ys_view = @view ys[inds]
        Ninvdiag_view = @view Ninvdiag[inds]

        @inbounds for (j, wtoa) in enumerate(wtoas)
            cwtoa = correct_toa(model, wtoa, params)
            dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)

            ys_view[j] = dphase / doppler_shifted_spin_frequency(cwtoa.toa_correction)
            ys_view[ntoas+j] = dm_residual(wtoa.dminfo, cwtoa.dm_correction)

            Ninvdiag_view[j] = 1.0 / scaled_toa_error_sqr(wtoa.toa, cwtoa.toa_correction)
            Ninvdiag_view[ntoas+j] =
                1.0 / scaled_dm_error_sqr(wtoa.dminfo, cwtoa.dm_correction)
        end
    end

    return ys, Ninvdiag
end
