using GeometricUnits
using Quadmath

export TOA

struct TOA
    val::GQ{Float128}
    err::GQ{Float64}
    freq::GQ{Float64}
    clkcorr::GQ{Float64}
    deltapn::GQ{Float64}
    ssb_obs_pos::Tuple{GQ{Float64}, GQ{Float64}, GQ{Float64}}
    ssb_obs_vel::Tuple{GQ{Float64}, GQ{Float64}, GQ{Float64}}
    obs_sun_pos::Tuple{GQ{Float64}, GQ{Float64}, GQ{Float64}}
    level::UInt
    phase::GQ{Float128}
end
