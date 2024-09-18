@testset "DispersionTaylor" begin
    toa = default_toa()
    ctoa = TOACorrection()

    tzrtoa = default_tzrtoa()
    ctzrtoa = TOACorrection()

    wtoa = default_wbtoa()
    cwtoa = WidebandTOACorrection()

    params = (DMEPOCH = time(53470.0 * day_to_s), DM = (GQ{-1}(4e16), GQ{-2}(1e11)))

    dmt = DispersionTaylor()
    @test dispersion_slope(dmt, toa, ctoa, params) == params.DM[1]
    @test delay(dmt, toa, ctoa, params) ==
          dispersion_slope(dmt, toa, ctoa, params) / toa.observing_frequency^Val(2)

    @test @ballocated(delay($dmt, $toa, $ctoa, $params)) == 0

    cwtoa1 = correct_toa(dmt, wtoa, cwtoa, params)
    @test cwtoa1.dm_correction.model_dm == dispersion_slope(dmt, toa, ctoa, params)

    display(dmt)
end

@testset "DispersionOffset" begin
    toa = default_toa()
    ctoa = TOACorrection()

    tzrtoa = default_tzrtoa()
    ctzrtoa = TOACorrection()

    wtoa = default_wbtoa()
    cwtoa = WidebandTOACorrection()

    params = (FDJUMPDM = (GQ{-1}(1e11), GQ{-1}(-1e11)),)

    jump_mask = BitMatrix([1 0 0; 0 1 0])
    dmoff = DispersionOffset(jump_mask)

    @test dispersion_slope(dmoff, tzrtoa, ctzrtoa, params) == GQ{-1}(0.0)
    @test dispersion_slope(dmoff, toa, ctoa, params) == -params.FDJUMPDM[1]

    @test @ballocated(delay($dmoff, $toa, $ctoa, $params)) == 0

    display(dmoff)

    jump_mask_ex = [1, 2, 0]
    dmoff_ex = ExclusiveDispersionOffset(jump_mask_ex)

    @test dispersion_slope(dmoff_ex, tzrtoa, ctzrtoa, params) == GQ{-1}(0.0)
    @test dispersion_slope(dmoff_ex, toa, ctoa, params) == -params.FDJUMPDM[1]

    toa1 = TOA(
        time(Double64(53470.0 * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        3,
    )
    ctoa1 = TOACorrection()
    @test dispersion_slope(dmoff_ex, toa1, ctoa1, params) ≈ GQ{-1}(0.0)

    display(dmoff_ex)
end

@testset "DispersionJump" begin
    toa = default_toa()
    ctoa = TOACorrection()

    tzrtoa = default_tzrtoa()
    ctzrtoa = TOACorrection()

    wtoa = default_wbtoa()
    cwtoa = WidebandTOACorrection()

    params = (DMJUMP = (GQ{-1}(1e11), GQ{-1}(-1e11)),)

    jump_mask = BitMatrix([1 0 0; 0 1 0])
    dmjump = DispersionJump(jump_mask)

    @test dispersion_slope(dmjump, tzrtoa, ctzrtoa, params) == GQ{-1}(0.0)
    @test dispersion_slope(dmjump, toa, ctoa, params) == -params.DMJUMP[1]

    @test delay(dmjump, toa, ctoa, params) == time(0.0)

    cwtoa1 = correct_toa(dmjump, wtoa, cwtoa, params)
    @test cwtoa1.toa_correction.delay == time(0.0)

    @test @ballocated(correct_toa($dmjump, $wtoa, $cwtoa, $params)) == 0

    display(dmjump)

    jump_mask_ex = [1, 2, 0]
    dmjump_ex = ExclusiveDispersionJump(jump_mask_ex)

    @test dispersion_slope(dmjump_ex, tzrtoa, ctzrtoa, params) == GQ{-1}(0.0)
    @test dispersion_slope(dmjump_ex, toa, ctoa, params) == -params.DMJUMP[1]

    toa1 = TOA(
        time(Double64(53470.0 * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        3,
    )
    ctoa1 = TOACorrection()
    @test dispersion_slope(dmjump_ex, toa1, ctoa1, params) ≈ GQ{-1}(0.0)

    @test delay(dmjump_ex, toa1, ctoa1, params) == time(0.0)

    display(dmjump_ex)
end
