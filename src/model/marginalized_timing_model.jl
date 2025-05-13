export MarginalizedTimingModel

"""The linearized and marginalized part of the timing model.
This represents small deviations of timing model parameters from their
default values when their effect on the residuals is (approximately) 
linear. Supported parameters include PHOFF, F, JUMP, and DMJUMP. This 
component doesn't have a `correct_toa()` method. The corresponding design
matrix is included in the `WoodburyKernel`.

Reference:
    [van Haasteren & Levin 2013](https://doi.org/10.1093/mnras/sts097)
"""
struct MarginalizedTimingModel <: Component
    prior_weights_inv::Vector{Float64}
    marginalized_param_names::Vector{String}

    function MarginalizedTimingModel(weights, pnames)
        @assert length(weights) == length(pnames)
        return new(1 ./ weights, pnames)
    end
end

is_gp_noise(::MarginalizedTimingModel) = true # COV_EXCL_LINE
get_gp_npars(mtm::MarginalizedTimingModel) = length(mtm.prior_weights_inv)
get_marginalized_param_names(mtm::MarginalizedTimingModel) = mtm.marginalized_param_names
calc_noise_weights_inv(mtm::MarginalizedTimingModel, ::NamedTuple) = mtm.prior_weights_inv
