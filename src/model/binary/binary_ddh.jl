export BinaryDDH

"""The Damour & Deruelle model for eccentric binaries with orthometric parametrization of
the Shapiro delay. Used for low to moderate inclination binaries.

References:
    [Damour & Deruelle 1986](https://ui.adsabs.harvard.edu/abs/1986AIHPA..44..263D/abstract), 
    [Freire & Wex 2010](http://doi.org/10.1111/j.1365-2966.2010.17319.x),
    [Weisberg & Huang 2016](http://doi.org/10.3847/0004-637X/829/1/55)
"""
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
