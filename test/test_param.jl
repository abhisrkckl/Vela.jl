@testset "parameter & param handler" begin
    pepoch = Parameter(:PEPOCH, time(56000.0 * day_to_s), true, "day", float(day_to_s))
    @test pepoch.name == :PEPOCH

    f0 = Parameter(:F0, frequency(100.0), false, "Hz", 1.0)
    f1 = Parameter(:F1, GQ{-2}(-1e-14), false, "Hz^2", 1.0)
    mparF = MultiParameter(:F, [f0, f1])
    @test length(mparF.parameters) == 2

    phoff = Parameter(:PHOFF, dimensionless(0.0), false, "", 1.0)

    ph = ParamHandler([pepoch, phoff], [mparF])
    @test get_free_param_names(ph) == ["PHOFF", "F0", "F1"]
    @test Set(keys(ph._default_params_tuple)) == Set([:PEPOCH, :PHOFF, :F])

    params = read_params(ph, [0.01, 100.01, -1.01e-14])
    @test Set(keys(params)) == Set([:PEPOCH, :PHOFF, :F])
    @test params.F == (frequency(100.01), GQ{-2}(-1.01e-14))

    @test all(get_scale_factors(ph) .> 0)

    @test length(read_param_values_to_vector(ph)) == length(get_free_param_names(ph))
    @test read_params(ph, read_param_values_to_vector(ph)) == ph._default_params_tuple
end
