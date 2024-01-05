using Vela
using Test
using GeometricUnits
using LinearAlgebra
using Quadmath
using JuliaFormatter

@testset "Vela" begin
    @testset "toa" begin
        ssb_obs_pos = distance.([18.0354099, 450.01472245, 195.05827732])
        ssb_obs_vel = speed.([-9.96231954e-05, 3.31555854e-06, 1.12968547e-06])
        obs_sun_pos = distance.([-15.89483533, -450.17185232, -195.18212616])

        ephem_vecs = EphemerisVectors(ssb_obs_pos, ssb_obs_vel, obs_sun_pos)

        # Aphelion distance
        @test sqrt(dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos)) < distance(509.0)
        @test sqrt(dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos)) < distance(509.0)

        # Angle
        @test acos(
            -dot(ephem_vecs.ssb_obs_pos, ephem_vecs.obs_sun_pos) / sqrt(
                dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos) *
                dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos),
            ),
        ) < 0.01

        # Speed of light
        @test dot(ephem_vecs.ssb_obs_vel, ephem_vecs.ssb_obs_vel) < 1

        toaval = time(parse(Float128, "4610197611.8464445127"))
        toaerr = time(1e-6)
        freq = frequency(1.4e9)
        phase = dimensionless(Float128(1000.0))

        toa1 = TOA(toaval, toaerr, freq, phase)
        @test is_barycentered(toa1)
        @test !toa1.tzr
        @test toa1.level == 0

        toa2 = TOA(toaval, toaerr, freq, phase, ephem_vecs)
        @test !is_barycentered(toa2)
        @test !toa2.tzr
        @test toa2.level == 0

        dt = time(1.0)
        toa3 = correct_toa_delay(toa1, dt)
        @test toa3.value == toa1.value - dt
        @test toa3.error == toa1.error
        @test toa3.frequency == toa1.frequency
        @test toa3.phase == toa1.phase
        @test is_barycentered(toa3) == is_barycentered(toa1)
        @test toa3.tzr == toa1.tzr
        @test toa3.level == 1

        dphi = dimensionless(0.3)
        toa4 = correct_toa_phase(toa3, dphi)
        @test toa4.value == toa3.value
        @test toa4.error == toa3.error
        @test toa4.frequency == toa3.frequency
        @test toa4.phase == toa3.phase + dphi
        @test is_barycentered(toa4) == is_barycentered(toa3)
        @test toa4.tzr == toa3.tzr
        @test toa4.level == 2
    end

    @testset "formatting" begin
        @test format(".")
    end
end
