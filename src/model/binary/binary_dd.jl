export BinaryDD

struct BinaryDD <: BinaryDDBase
    use_fbx::Bool
end

shapiro_delay_params(::BinaryDD, params::NamedTuple) = params.M2, params.SINI

function show(io::IO, dd::BinaryDD)
    mode = dd.use_fbx ? "FBX" : "PB"
    print(io, "BinaryDD(with $mode)")
end
