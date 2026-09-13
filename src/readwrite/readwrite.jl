function save_pulsar_data(
    filename::String,
    model::TimingModel,
    toas::Vector{T},
) where {T<:TOABase}
    JLSO.save(filename, :model => model, :toas => toas)
end

function load_pulsar_data(filename::String)
    data = JLSO.load(filename)
    model, toas = data[:model], data[:toas]

    if isa(model.kernel, WoodburyKernel)
        # The WoodburyKernel must be re-made because the kernel workspace size
        # depends on the number of available threads.
        kernel = WoodburyKernel(
            model.kernel.inner_kernel,
            model.kernel.gp_components,
            model.kernel.noise_basis,
        )
        return TimingModel(
            model.pulsar_name,
            model.ephem,
            model.clock,
            model.units,
            model.epoch,
            model.components,
            kernel,
            model.param_handler,
            model.tzr_toa,
            model.priors,
        ),
        toas
    else
        return model, toas
    end
end
