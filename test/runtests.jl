using Vela
using Test
using GeometricUnits
using LinearAlgebra
using DoubleFloats
using JuliaFormatter
using BenchmarkTools
using JLSO
using Distributions

const day_to_s = 86400

const ssb_obs_pos = distance.((18.0354099, 450.01472245, 195.05827732))
const ssb_obs_vel = speed.((-9.96231954e-05, 3.31555854e-06, 1.12968547e-06))
const obs_sun_pos = distance.((-15.89483533, -450.17185232, -195.18212616))
const obs_jupiter_pos = distance.((-1610.1796849, -1706.87348483, -681.22381513))
const obs_saturn_pos = distance.((-2392.85651431, 3109.13626083, 1405.71912274))
const obs_venus_pos = distance.((140.85922773, -217.65571843, -74.64804201))
const obs_uranus_pos = distance.((9936.62957939, -3089.07377113, -1486.17339104))
const obs_neptune_pos = distance.((11518.60924426, -9405.0693235, -4126.91030657))
const obs_earth_pos = distance.((0.01199435, 0.01159591, -0.01316261))

@testset "Vela" verbose = true begin

    include("test_ephemeris.jl")

    include("test_toa.jl")

    include("test_param.jl")

    # include("test_components.jl")

    # include("test_priors.jl")

    # include("test_pure_rotator.jl")

    # include("test_NGC6440E.jl")

    @testset "formatting" begin
        @test format(Vela)
    end
end
