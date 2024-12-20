@testset "ChromaticTaylor" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        TNCHROMIDX = dimensionless(2.0),
        CMEPOCH = time((53470.0 - epoch_mjd) * day_to_s),
        CM = (GQ{1}(4e4), GQ{0}(1e-1)),
        DMEPOCH = time((53470.0 - epoch_mjd) * day_to_s),
        DM = (GQ{-1}(4e16), GQ{-2}(1e11)),
    )

    cmt = ChromaticTaylor()
    dmt = DispersionTaylor()
    @test chromatic_slope(cmt, toa, ctoa, params) == params.CM[1]

    # When TNCHROMIDX == 2, DM and CM that have equal numerical values
    # in the par file must give equal delays. In our units, this corresponds
    # to `value(DM) / value(CM) == 1e12`. See the `params` tuple above.
    @test delay(dmt, toa, ctoa, params) == delay(cmt, toa, ctoa, params)

    @test @ballocated(delay($cmt, $toa, $ctoa, $params)) == 0
end

@testset "ChromaticPiecewise" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (TNCHROMIDX = dimensionless(4.0), CMX_ = (GQ{-1}(4e12), GQ{-1}(1e12)))

    mask = UInt[0]
    cmx = ChromaticPiecewise(mask)
    @test chromatic_slope(cmx, toa, ctoa, params) == GQ{-1}(0.0)
    @test @ballocated(delay($cmx, $toa, $ctoa, $params)) == 0

    mask = UInt[1]
    cmx = ChromaticPiecewise(mask)
    @test chromatic_slope(cmx, toa, ctoa, params) == params.CMX_[1]
    @test @ballocated(delay($cmx, $toa, $ctoa, $params)) == 0

    tzrtoa = default_tzrtoa()
    ctzrtoa = TOACorrection()
    @test chromatic_slope(cmx, tzrtoa, ctzrtoa, params) == GQ{-1}(0.0)
end
