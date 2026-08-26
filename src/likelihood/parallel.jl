function calc_noise_weights_inv(
    kernel::WoodburyKernel,
    params_list::AbstractVector{<:NamedTuple},
)
    return map(params -> calc_noise_weights_inv(kernel, params), params_list)
end

function calc_resids_and_Ninvdiag(
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

function calc_resids_and_Ninvdiag(
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
