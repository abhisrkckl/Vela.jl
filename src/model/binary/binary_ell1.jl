export BinaryELL1

"""A binary model representing a nearly circular orbit."""
struct BinaryELL1 <: BinaryELL1Base
    use_fbx::Bool
end

shapiro_delay_params(::BinaryELL1, params::NamedTuple) = params.M2, params.SINI
