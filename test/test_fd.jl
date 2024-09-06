@testset "FrequencyDependent" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    params = (FD = (time(0.1), time(0.2)),)

    fd = FrequencyDependent()
    @test isfinite(delay(fd, ctoa, params))
    @test @ballocated(delay($fd, $ctoa, $params)) == 0
end
