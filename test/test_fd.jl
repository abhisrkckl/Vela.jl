@testset "FrequencyDependent" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    params = (FD = (time(0.1), time(0.2)),)

    fd = FrequencyDependent()
    @test isfinite(delay(fd, ctoa, params))
    @test @ballocated(delay($fd, $ctoa, $params)) == 0
end

@testset "FrequencyDependentJump" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    params = (FDJUMP = (time(0.1), time(0.2), time(0.015), time(0.031)),)

    mask = BitMatrix([
        1 1 0 0
        0 0 1 1
        1 1 0 0
        0 0 1 1
    ])
    exponents = UInt[1, 2, 1, 2]
    fdjump = FrequencyDependentJump(mask, exponents)

    @test delay(fdjump, ctoa, params) > time(0.0)
    @test @ballocated(delay($fdjump, $ctoa, $params)) == 0

    tzrtoa = default_tzrtoa()
    ctzrtoa = CorrectedTOA(tzrtoa)
    @test delay(fdjump, ctzrtoa, params) == time(0.0)
    @test @ballocated(delay($fdjump, $ctzrtoa, $params)) == 0

    display(fdjump)
end
