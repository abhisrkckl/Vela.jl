@testset "parameter & param handler" begin
    pepoch = Parameter(:PEPOCH, time(56000.0 * day_to_s), true, "day", day_to_s)
    @test pepoch.name == :PEPOCH
    @test_throws AssertionError read_param(pepoch, 56000.0 * day_to_s)

    f0 = Parameter(:F0, frequency(100.0), false, "Hz", 1.0)
    @test read_param(f0, 101.0) == frequency(101.0)

    f1 = Parameter(:F1, GQ(-1e-14, -2), false, "Hz^2", 1.0)
    @test read_param(f1, -1.1e-14) == GQ(-1.1e-14, -2)

    mparF = MultiParameter(:F, [f0, f1])
    @test length(mparF.parameters) == 2

    phoff = Parameter(:PHOFF, dimensionless(0.0), false, "", 1.0)

    ph = ParamHandler([pepoch, phoff], [mparF])
    @test get_free_param_names(ph) == ["PHOFF", "F0", "F1"]
    @test Set(keys(ph._default_params_tuple)) == Set([:PEPOCH, :PHOFF, :F])

    params = read_params(ph, [0.01, 100.01, -1.01e-14])
    @test Set(keys(params)) == Set([:PEPOCH, :PHOFF, :F])
    @test params.F == (frequency(100.01), GQ(-1.01e-14, -2))

    @test all(get_scale_factors(ph) .> 0)
end
