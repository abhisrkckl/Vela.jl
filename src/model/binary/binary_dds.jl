export BinaryDDS

"""The Damour & Deruelle model for eccentric binaries with an alternative parametrization of
Shapiro delay applicable to almost edge-on orbits.

References:
    [Damour & Deruelle 1986](https://ui.adsabs.harvard.edu/abs/1986AIHPA..44..263D/abstract),
    [Kramer+ 2006](http://doi.org/10.1126/science.1132305),
    [Rafikov & Lai 2006](http://doi.org/10.1103/PhysRevD.73.063003)
"""
struct BinaryDDS <: BinaryDDBase
    use_fbx::Bool
end

function shapiro_delay_params(::BinaryDDS, params::NamedTuple)
    smax = params.SHAPMAX
    m2 = params.M2
    sini = 1 - exp(-smax)
    return m2, sini
end
