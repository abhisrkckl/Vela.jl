@testset "toa" begin
    toaval = time(parse(Double64, "4610197611.8464445127"))
    toaerr = time(1e-6)
    freq = frequency(1.4e9)
    pulse_number = dimensionless(Double64(1000.0))

    ephem = SolarSystemEphemeris(
        ssb_obs_pos,
        ssb_obs_vel,
        obs_sun_pos,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
        obs_earth_pos,
    )

    # TOA value should be of type GQ{Double64}.
    @test_throws MethodError TOA(time(4610197611.8), toaerr, freq, pulse_number, ephem, 1)

    # Wrong dimensions for TOA value.
    @test_throws MethodError TOA(
        dimensionless(toaval.x),
        toaerr,
        freq,
        pulse_number,
        ephem,
        1,
    )

    # Wrong dimensions for TOA error.
    @test_throws MethodError TOA(toaval, dimensionless(1e-6), freq, pulse_number, ephem, 1)

    # Wrong dimensions for TOA observing_frequency.
    @test_throws MethodError TOA(toaval, toaerr, time(1.4e9), phase, ephem, 1)

    # Wrong dimensions for TOA pulse_number.
    @test_throws MethodError TOA(toaval, toaerr, freq, time(1000.0), ephem, 1)

    toa1 = TOA(toaval, toaerr, freq, pulse_number, ephem, 1)
    @test !is_tzr(toa1)

    toacorr1 = TOACorrection()

    dt = time(1.0)
    toacorr2 = correct_toa_delay(toacorr1; delay = dt)
    @test toacorr2.delay == toacorr1.delay + dt
    @test corrected_toa_value(toa1, toacorr2) == corrected_toa_value(toa1, toacorr1) - dt
    @test toacorr2.phase == toacorr1.phase
    @test toacorr2.efac == toacorr1.efac
    @test toacorr2.equad2 == toacorr1.equad2
    @test toacorr2.doppler == toacorr1.doppler
    @test is_barycentered(toacorr2) == is_barycentered(toacorr1)

    dphi = dimensionless(0.3)
    toacorr3 = correct_toa_phase(toacorr2; phase = dphi)
    @test toacorr3.delay == toacorr2.delay
    @test toacorr3.phase == toacorr2.phase + dphi
    @test phase_residual(toa1, toacorr3) == phase_residual(toa1, toacorr2) + dphi
    @test toacorr3.efac == toacorr2.efac
    @test toacorr3.equad2 == toacorr2.equad2
    @test toacorr3.doppler == toacorr2.doppler
    @test is_barycentered(toacorr3) == is_barycentered(toacorr2)

    efac = dimensionless(1.1)
    equad2 = time(1e-6)^Val(2)
    toacorr4 = correct_toa_error(toacorr3; efac = efac, equad2 = equad2)
    @test toacorr4.delay == toacorr3.delay
    @test toacorr4.phase == toacorr3.phase
    @test toacorr4.efac == toacorr3.efac * efac
    @test toacorr4.equad2 == toacorr3.equad2 + equad2
    @test scaled_toa_error_sqr(toa1, toacorr4) ≈
          (scaled_toa_error_sqr(toa1, toacorr3) + equad2) * efac^2
    @test toacorr4.doppler == toacorr3.doppler
    @test is_barycentered(toacorr4) == is_barycentered(toacorr3)

    @testset "tzr_toa" begin
        tzrtoa = make_tzr_toa(toaval, freq, ephem)
        @test is_tzr(tzrtoa)
        @test !is_barycentered(tzrtoa)
        @test tzrtoa.error == time(0.0)
        @test tzrtoa.index == 0
    end
end
