export BinaryELL1

"""A binary model representing a nearly circular orbit.

Reference:
    [Lange+ 2001](http://doi.org/10.1046/j.1365-8711.2001.04606.x)
"""
struct BinaryELL1 <: BinaryELL1Base
    use_fbx::Bool
end

shapiro_delay_params(::BinaryELL1, params::NamedTuple) = params.M2, params.SINI
