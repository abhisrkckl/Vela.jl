export BinaryDDH

struct BinaryDDH <: BinaryDDBase
    use_fbx::Bool
end

function shapiro_delay_params(::BinaryDDH, params::NamedTuple)
    h3 = params.H3
    Ϛ = params.STIGMA
    m2 = h3 / Ϛ^3
    sini = 2 * Ϛ / (1 + Ϛ^2)
    return m2, sini
end
