export BinaryDD

"""
    BinaryDD
    
The Damour & Deruelle binary model for eccentric binaries.

Reference:
    [Damour & Deruelle 1986](https://ui.adsabs.harvard.edu/abs/1986AIHPA..44..263D/abstract)
"""
struct BinaryDD <: BinaryDDBase
    use_fbx::Bool
end

shapiro_delay_params(::BinaryDD, params::NamedTuple) = params.M2, params.SINI

function show(io::IO, dd::BinaryDD)
    mode = dd.use_fbx ? "FBX" : "PB"
    print(io, "BinaryDD(with $mode)")
end
