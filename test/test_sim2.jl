@testset "sim2" begin
    model, toas = Vela.load_pulsar_data("datafiles/sim2.jlso")

    @testset "file save" begin
        Vela.save_pulsar_data("__test_sim2.jlso", model, toas)
        @test isfile("__test_sim2.jlso")
    end

    @testset "copy" begin
        m2 = TimingModel(
            model.pulsar_name,
            model.ephem,
            model.clock,
            model.units,
            model.components,
            model.kernel,
            model.param_handler,
            model.tzr_toa,
            model.priors,
        )
    end

    @testset "read_toas" begin
        @test !any([is_tzr(toa) for toa in toas])
        @test length(toas) == 2000
        @test all([
            frequency(4.9e8) < toa.observing_frequency < frequency(1.51e9) for toa in toas
        ])
        @test all([
            time(52999.0 * day_to_s) < toa.value < time(57001.0 * day_to_s) for toa in toas
        ])
        @test all([modf(toa.pulse_number.x)[1] == 0 for toa in toas])
        @test all([toa.error > time(0.0) for toa in toas])
    end

    @testset "tzr_toa" begin
        tzrtoa = model.tzr_toa
        @test is_tzr(tzrtoa)
        @test tzrtoa.error > time(0.0)
        @test tzrtoa.pulse_number == dimensionless(0.0)
        @test frequency(4.9e8) < tzrtoa.observing_frequency < frequency(1.51e9)
        @test time(52999.0 * day_to_s) < tzrtoa.value < time(57001.0 * day_to_s)
    end

    @testset "param_handler" begin
        param_handler = model.param_handler
        @test Set(get_free_param_names(param_handler)) ==
              Set(["F0", "F1", "PHOFF", "RAJ", "DECJ", "EFAC1", "ECORR1"])
        @test length(param_handler.multi_params) + length(param_handler.single_params) ==
              length(param_handler._default_params_tuple)
        @test length(get_free_param_names(param_handler)) ==
              length(param_handler._free_indices)
        @test sizeof(param_handler._default_params_tuple) ==
              sizeof(GQ{0,Float64}) * length(param_handler._default_values)

        @test length(
            read_param_values_to_vector(
                model.param_handler,
                model.param_handler._default_params_tuple,
            ),
        ) == length(param_handler._free_indices)

        @test length(
            read_param_values_to_vector(model, model.param_handler._default_params_tuple),
        ) == length(param_handler._free_indices)

        @test all(isfinite.(get_scale_factors(model))) &&
              !any(iszero.(get_scale_factors(model)))
    end

    @testset "components" begin
        components = model.components
        @test length(components) == 5

        @test isa(model.kernel, EcorrKernel)

        @test isa(components[1], SolarSystem)
        @test !components[1].ecliptic_coordinates
        @test !components[1].planet_shapiro

        @test isa(components[2], DispersionTaylor)

        @test isa(components[3], Spindown)

        @test isa(components[4], PhaseOffset)

        @test isa(components[5], MeasurementNoise)
    end

    @testset "form_residuals" begin
        params = model.param_handler._default_params_tuple
        res = form_residuals(model, toas, params)
        @test all(abs(r) < 7 * toa.error for (r, toa) in zip(res, toas))
    end

    @testset "calc_chi2" begin
        calc_chi2_s = get_chi2_serial_func(model, toas)
        calc_chi2_p = get_chi2_parallel_func(model, toas)
        parv = read_param_values_to_vector(model.param_handler)
        # parnp = PyArray(parv)
        chi2_s = calc_chi2_s(parv)
        chi2_p = calc_chi2_p(parv)
        @test chi2_s / degrees_of_freedom(model, toas) < 1.2
        @test chi2_s ≈ chi2_p
    end

    @testset "calc_lnlike" begin
        calc_lnlike_s = get_lnlike_serial_func(model, toas)
        calc_lnlike_p = get_lnlike_parallel_func(model, toas)
        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model)
        # parnp = PyArray(parv)
        @test calc_lnlike_s(parv) ≈ calc_lnlike_p(parv)

        @test @ballocated(Vela.calc_lnlike_serial($model, $toas, $params)) == 0

        parv1 = read_param_values_to_vector(model.param_handler, params)
        parv1[end] *= 2
        @test calc_lnlike(model, toas, parv1) < calc_lnlike(model, toas, params)
    end

    @testset "priors" begin
        calc_lnprior = get_lnprior_func(model)
        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model.param_handler, params)
        @test isfinite(calc_lnprior(model.param_handler._default_params_tuple))
        @test calc_lnprior(params) == calc_lnprior(parv)

        calc_lnpost = get_lnpost_func(model, toas)
        @test isfinite(calc_lnpost(params))

        prior_transform = get_prior_transform_func(model)
        halfs = fill(0.5, length(parv))
        @test all(isfinite.(prior_transform(halfs)))
        # @test all(prior_transform(halfs) .≈ parv)
    end
end
