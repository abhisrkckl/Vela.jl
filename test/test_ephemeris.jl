@testset "ephemeris" begin
    # Wrong dimensions
    @test_throws MethodError SolarSystemEphemeris(
        ssb_obs_vel,
        ssb_obs_vel,
        obs_sun_pos,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
        # obs_earth_pos,
    )

    # Wrong dimensions
    @test_throws MethodError SolarSystemEphemeris(
        ssb_obs_pos,
        ssb_obs_pos,
        obs_sun_pos,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
        # obs_earth_pos,
    )

    # Wrong dimensions
    @test_throws MethodError SolarSystemEphemeris(
        ssb_obs_pos,
        ssb_obs_vel,
        ssb_obs_vel,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
        # obs_earth_pos,
    )

    # ssb_obs_vel is too large.
    @test_throws AssertionError SolarSystemEphemeris(
        ssb_obs_pos,
        1e6 .* ssb_obs_vel,
        obs_sun_pos,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
        # obs_earth_pos,
    )

    ephem_vecs = SolarSystemEphemeris(
        ssb_obs_pos,
        ssb_obs_vel,
        obs_sun_pos,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
        # obs_earth_pos,
    )

    # ssb_obs_pos and obs_sun_pos should be less than the Aphelion distance
    @test sqrt(dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos)) < distance(509.0)
    @test sqrt(dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos)) < distance(509.0)

    # ssb_obs_pos and obs_sun_pos should have a small angle.
    @test acos(
        -dot(ephem_vecs.ssb_obs_pos, ephem_vecs.obs_sun_pos) / sqrt(
            dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos) *
            dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos),
        ),
    ) < 0.01

    # |ssb_obs_vel| should be less than the speed of light
    @test dot(ephem_vecs.ssb_obs_vel, ephem_vecs.ssb_obs_vel) < 1

    @test sizeof(ephem_vecs) == 8 * 3 * sizeof(GQ{0,Float64})
end
