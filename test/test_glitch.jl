@testset "Glitch" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    toa2 = TOA(
        time(Double64(54954.0 * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        false,
        default_ephem(),
        1,
    )
    ctoa2 = CorrectedTOA(toa2)

    params = (
        GLEP_ = (time(day_to_s * 54952.92239),),
        GLPH_ = (dimensionless(0.3989),),
        GLF0_ = (frequency(1.750 - 06),),
        GLF1_ = (GQ{-2}(-6.572e-15),),
        GLF2_ = (GQ{-3}(0.0),),
        GLF0D_ = (frequency(0.0),),
        GLTD_ = (time(0.0),),
    )

    glitch = Glitch()

    @test phase(glitch, ctoa, params) == dimensionless(0.0)
    @test phase(glitch, ctoa2, params) != dimensionless(0.0)

    @test @ballocated(phase($glitch, $ctoa2, $params)) == 0
end
