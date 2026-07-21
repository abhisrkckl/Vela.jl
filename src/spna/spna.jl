export make_SPNA

abstract type SPNABase end

struct NarrowbandSPNA{
    KernelType<:WoodburyKernel,
    ComponentsTuple<:Tuple,
    PriorsTuple<:Tuple,
} <: SPNABase
    kernel::KernelType
    components::ComponentsTuple
    resids::Vector{NarrowbandResid}
    param_handler::ParamHandler
    priors::PriorsTuple
end

NarrowbandSPNA(model::TimingModel, toas::Vector{TOA}, Mtm::Array{Float64,2}) =
    NarrowbandSPNA(
        make_spna_kernel(model, toas, Mtm),
        filter(c -> isa(c, WhiteNoiseComponent), model.components),
        make_resids(model, toas),
        make_spna_param_handler(model),
        model.priors[(get_num_timing_params(model)+1):end],
    )

struct WidebandSPNA{KernelType<:WoodburyKernel,ComponentsTuple<:Tuple,PriorsTuple<:Tuple} <:
       SPNABase
    kernel::KernelType
    components::ComponentsTuple
    resids::Vector{WidebandResid}
    param_handler::ParamHandler
    priors::PriorsTuple
end

WidebandSPNA(model::TimingModel, wtoas::Vector{WidebandTOA}, Mtm::Array{Float64,2}) =
    WidebandSPNA(
        make_spna_kernel(model, wtoas, Mtm),
        filter(c -> isa(c, WhiteNoiseComponent), model.components),
        make_resids(model, wtoas),
        make_spna_param_handler(model),
        model.priors[(get_num_timing_params(model)+1):end],
    )

make_SPNA(model::TimingModel, toas::Vector{TOA}, Mtm::Array{Float64,2}) =
    NarrowbandSPNA(model, toas, Mtm)
make_SPNA(model::TimingModel, wtoas::Vector{WidebandTOA}, Mtm::Array{Float64,2}) =
    WidebandSPNA(model, wtoas, Mtm)

show(io::IO, spna::NarrowbandSPNA) =
    print(io, "NarrowbandSPNA[$(size(spna.kernel.noise_basis))]")
show(io::IO, spna::WidebandSPNA) =
    print(io, "WidebandSPNA[$(size(spna.kernel.noise_basis))]")
show(io::IO, ::MIME"text/plain", spna::SPNABase) = show(io, spna)
