function make_spna_kernel(model::TimingModel, toas::Vector{TOA}, Mtm::Array{Float64,2})
    ntoas, ntmpars = size(Mtm)
    @assert ntoas == length(toas) "Number of TOAs must match number of rows in the design matrix."
    @assert ntmpars == get_num_timing_params(model) "Number of timing model parameters must match number of columns in the design matrix."

    timing_param_names = get_free_param_names(model)[1:ntmpars]
    prior_weights_inv = fill(1e-40, ntmpars)
    mtm = MarginalizedTimingModel(prior_weights_inv, timing_param_names)

    if isa(model.kernel, WoodburyKernel)
        inner_kernel = model.kernel.inner_kernel
        M = hcat(Mtm, model.kernel.noise_basis)
        gp_components = (mtm, model.kernel.gp_components...)
    else
        inner_kernel = model.kernel
        M = Mtm
        gp_components = (mtm,)
    end

    return WoodburyKernel(inner_kernel, gp_components, M)
end

function make_spna_param_handler(model::TimingModel)::ParamHandler
    single_params = filter(p -> p.noise, model.param_handler.single_params)
    multi_params = filter(mp -> mp.parameters[1].noise, model.param_handler.multi_params)
    return ParamHandler(single_params, multi_params)
end
