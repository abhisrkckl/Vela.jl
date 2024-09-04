export BinaryDDS

struct BinaryDDS <: BinaryDDBase
    use_fbx::Bool
end

function shapiro_delay_params(::BinaryDDS, params::NamedTuple)
    smax = params.SHAPMAX
    m2 = params.M2
    sini = 1 - exp(-smax)
    return m2, sini
end
