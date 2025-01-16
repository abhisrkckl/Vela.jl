struct WoodburyKernel{InnerKernel<:Kernel,SectionsTuple<:Tuple} <: Kernel
    inner_kernel::InnerKernel
    sections::SectionsTuple
    basis::Matrix{Float64}
end

abstract type KernelSection end
