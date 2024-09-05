@testset "wideband toa" begin
    toaval = time(parse(Double64, "4610197611.8464445127"))
    toaerr = time(1e-6)
    freq = frequency(1.4e9)
    pulse_number = dimensionless(Double64(1000.0))
    barycentered = false

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

    toa = TOA(toaval, toaerr, freq, pulse_number, barycentered, ephem, 1)
    dminfo = DMInfo(GQ{-1}(1e16), GQ{-1}(1e11))
    wtoa = WidebandTOA(toa, dminfo)
    cwtoa1 = CorrectedWidebandTOA(wtoa)

    @test cwtoa1.corrected_toa.level == 0

    @test dm_residual(cwtoa1.corrected_dminfo) == dminfo.value
    @test scaled_dm_error_sqr(cwtoa1.corrected_dminfo) == dminfo.error^Val(2)

    ddm = GQ{-1}(1e14)
    cdminfo = cwtoa1.corrected_dminfo
    cdminfo2 = correct_dminfo(
        cdminfo;
        delta_dm = ddm,
        dmefac = dimensionless(1.1),
        dmequad2 = GQ{-2}(1e-12),
    )
    @test cdminfo2.model_dm == cdminfo.model_dm + ddm
    @test cdminfo2.dmefac == dimensionless(1.1)
    @test cdminfo2.dmequad2 == GQ{-2}(1e-12)
    # @test corrected_toa_value(ctoa2) == corrected_toa_value(ctoa1) - dt
    # @test ctoa2.phase == ctoa1.phase
    # @test ctoa2.efac == ctoa1.efac
    # @test ctoa2.equad2 == ctoa1.equad2
    # @test ctoa2.doppler == ctoa1.doppler
    # @test ctoa2.barycentered == ctoa1.barycentered
    # @test ctoa2.level == 1

    # dphi = dimensionless(0.3)
    # ctoa3 = correct_toa(ctoa2; phase = dphi)
    # @test ctoa3.delay == ctoa2.delay
    # @test ctoa3.phase == ctoa2.phase + dphi
    # @test phase_residual(ctoa3) == phase_residual(ctoa2) + dphi
    # @test ctoa3.efac == ctoa2.efac
    # @test ctoa3.equad2 == ctoa2.equad2
    # @test ctoa3.doppler == ctoa2.doppler
    # @test ctoa3.barycentered == ctoa2.barycentered
    # @test ctoa3.level == 2

    # efac = dimensionless(1.1)
    # equad2 = time(1e-6)^Val(2)
    # ctoa4 = correct_toa(ctoa3; efac = efac, equad2 = equad2)
    # @test ctoa4.delay == ctoa3.delay
    # @test ctoa4.phase == ctoa3.phase
    # @test ctoa4.efac == ctoa3.efac * efac
    # @test ctoa4.equad2 == ctoa3.equad2 + equad2
    # @test scaled_toa_error_sqr(ctoa4) ≈ (scaled_toa_error_sqr(ctoa3) + equad2) * efac^2
    # @test ctoa4.doppler == ctoa3.doppler
    # @test ctoa4.barycentered == ctoa3.barycentered
    # @test ctoa4.level == 3

    # @testset "tzr_toa" begin
    #     tzrtoa = make_tzr_toa(toaval, freq, true, ephem)
    #     @test tzrtoa.tzr
    #     @test tzrtoa.barycentered
    #     @test tzrtoa.error == time(0.0)
    #     @test tzrtoa.index == 0

    #     ctzrtoa = CorrectedTOA(tzrtoa)
    #     @test ctzrtoa.level == 0
    # end
end