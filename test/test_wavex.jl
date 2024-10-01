@testset "WaveX" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        WXEPOCH = time((53470.0 - epoch_mjd) * day_to_s),
        WXFREQ_ = (frequency(1e-9), frequency(2e-9), frequency(3e-9)),
        WXSIN_ = (time(1.2e-6), time(5.1e-7), time(2.5e-7)),
        WXCOS_ = (time(-1.3e-6), time(5.2e-7), time(2.6e-7)),
    )

    wx = WaveX()
    @test isfinite(delay(wx, toa, ctoa, params))
    @test @ballocated(delay($wx, $toa, $ctoa, $params)) == 0
end

@testset "DMWaveX" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        DMWXEPOCH = time((53470.0 - epoch_mjd) * day_to_s),
        DMWXFREQ_ = (frequency(1e-9), frequency(2e-9), frequency(3e-9)),
        DMWXSIN_ = (GQ{-1}(1.2e+4), GQ{-1}(5.1e+3), GQ{-1}(2.5e+2)),
        DMWXCOS_ = (GQ{-1}(-1.3e+4), GQ{-1}(5.2e+3), GQ{-1}(2.6e+3)),
    )

    dmwx = DMWaveX()
    @test isfinite(delay(dmwx, toa, ctoa, params))
    @test @ballocated(delay($dmwx, $toa, $ctoa, $params)) == 0
end

@testset "CMWaveX" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        CMWXEPOCH = time((53470.0 - epoch_mjd) * day_to_s),
        CMWXFREQ_ = (frequency(1e-9), frequency(2e-9), frequency(3e-9)),
        CMWXSIN_ = (GQ{1}(1.2e+4), GQ{1}(5.1e+3), GQ{1}(2.5e+2)),
        CMWXCOS_ = (GQ{1}(-1.3e+4), GQ{1}(5.2e+3), GQ{1}(2.6e+3)),
        TNCHROMIDX = dimensionless(2.0),
    )

    cmwx = CMWaveX()
    @test isfinite(delay(cmwx, toa, ctoa, params))
    @test @ballocated(delay($cmwx, $toa, $ctoa, $params)) == 0
end
