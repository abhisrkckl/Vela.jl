@testset "DispersionTaylor" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    tzrtoa = default_tzrtoa()
    ctzrtoa = CorrectedTOA(tzrtoa)

    wtoa = default_wbtoa()
    cwtoa = CorrectedWidebandTOA(wtoa)

    params = (DMEPOCH = time(53470.0 * day_to_s), DM = (GQ{-1}(4e16), GQ{-2}(1e11)))

    dmt = DispersionTaylor()
    @test dispersion_slope(dmt, ctoa, params) == params.DM[1]
    @test delay(dmt, ctoa, params) ==
          dispersion_slope(dmt, ctoa, params) / ctoa.toa.observing_frequency^Val(2)

    @test @ballocated(delay($dmt, $ctoa, $params)) == 0

    cwtoa1 = correct_toa(dmt, cwtoa, params)
    @test cwtoa1.corrected_dminfo.model_dm == dispersion_slope(dmt, ctoa, params)

    display(dmt)
end

@testset "DispersionOffset" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    tzrtoa = default_tzrtoa()
    ctzrtoa = CorrectedTOA(tzrtoa)

    wtoa = default_wbtoa()
    cwtoa = CorrectedWidebandTOA(wtoa)

    params = (FDJUMPDM = (GQ{-1}(1e11), GQ{-1}(-1e11)),)

    jump_mask = BitMatrix([1 0 0; 0 1 0])
    dmoff = DispersionOffset(jump_mask)

    @test dispersion_slope(dmoff, ctzrtoa, params) == GQ{-1}(0.0)
    @test dispersion_slope(dmoff, ctoa, params) == -params.FDJUMPDM[1]

    @test @ballocated(delay($dmoff, $ctoa, $params)) == 0

    display(dmoff)

    jump_mask_ex = [1, 2, 0]
    dmoff_ex = ExclusiveDispersionOffset(jump_mask_ex)

    @test dispersion_slope(dmoff_ex, ctzrtoa, params) == GQ{-1}(0.0)
    @test dispersion_slope(dmoff_ex, ctoa, params) == -params.FDJUMPDM[1]

    toa1 = TOA(
        time(Double64(53470.0 * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        3,
    )
    ctoa1 = CorrectedTOA(toa1)
    @test dispersion_slope(dmoff_ex, ctoa1, params) ≈ GQ{-1}(0.0)

    display(dmoff_ex)
end

@testset "DispersionJump" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    tzrtoa = default_tzrtoa()
    ctzrtoa = CorrectedTOA(tzrtoa)

    wtoa = default_wbtoa()
    cwtoa = CorrectedWidebandTOA(wtoa)

    params = (DMJUMP = (GQ{-1}(1e11), GQ{-1}(-1e11)),)

    jump_mask = BitMatrix([1 0 0; 0 1 0])
    dmjump = DispersionJump(jump_mask)

    @test dispersion_slope(dmjump, ctzrtoa, params) == GQ{-1}(0.0)
    @test dispersion_slope(dmjump, ctoa, params) == -params.DMJUMP[1]

    @test delay(dmjump, ctoa, params) == time(0.0)

    cwtoa1 = correct_toa(dmjump, cwtoa, params)
    @test cwtoa1.corrected_toa.delay == time(0.0)

    @test @ballocated(correct_toa($dmjump, $cwtoa, $params)) == 0

    display(dmjump)

    jump_mask_ex = [1, 2, 0]
    dmjump_ex = ExclusiveDispersionJump(jump_mask_ex)

    @test dispersion_slope(dmjump_ex, ctzrtoa, params) == GQ{-1}(0.0)
    @test dispersion_slope(dmjump_ex, ctoa, params) == -params.DMJUMP[1]

    toa1 = TOA(
        time(Double64(53470.0 * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        3,
    )
    ctoa1 = CorrectedTOA(toa1)
    @test dispersion_slope(dmjump_ex, ctoa1, params) ≈ GQ{-1}(0.0)

    @test delay(dmjump_ex, ctoa1, params) == time(0.0)

    display(dmjump_ex)
end
