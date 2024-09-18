@testset "ChromaticTaylor" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        TNCHROMIDX = dimensionless(2.0),
        CMEPOCH = time(53470.0 * day_to_s),
        CM = (GQ{1}(4e4), GQ{0}(1e-1)),
        DMEPOCH = time(53470.0 * day_to_s),
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
