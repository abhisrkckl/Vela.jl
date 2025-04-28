@testset "sim_sw.wb" begin
    model, wtoas = Vela.load_pulsar_data("datafiles/sim_sw.wb.jlso")

    pepoch_mjd = value(model.epoch) / day_to_s

    @testset "file save" begin
        Vela.save_pulsar_data("__test_wb.jlso", model, wtoas)
        @test isfile("__test_wb.jlso")
    end

    @testset "repr" begin
        @test startswith(string(wtoas[1]), "WidebandTOA")
        display(wtoas)
        display(wtoas[1])
    end

    @testset "read_toas" begin
        @test !any([is_tzr(wtoa.toa) for wtoa in wtoas])
        @test length(wtoas) == 500
        @test all([
            frequency(1.3e9) < wtoa.toa.observing_frequency < frequency(1.5e9) for
            wtoa in wtoas
        ])
        @test all([
            time((53999.0 - pepoch_mjd) * day_to_s) <
            wtoa.toa.value <
            time((56001.0 - pepoch_mjd) * day_to_s) for wtoa in wtoas
        ])
        @test all([modf(wtoa.toa.pulse_number.x)[1] == 0 for wtoa in wtoas])
        @test all([wtoa.toa.error > time(0.0) for wtoa in wtoas])
    end

    @testset "param_handler" begin
        param_handler = model.param_handler
        @test Set(get_free_param_names(param_handler)) == Set([
            "F0",
            "F1",
            "PHOFF",
            "ELONG",
            "ELAT",
            "DM",
            "DM1",
            "NE_SW",
            "DMJUMP1",
            "DMEFAC1",
            "DMEQUAD1",
        ])
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
        @test length(components) == 7

        @test isa(components[1], SolarSystem)
        @test components[1].ecliptic_coordinates
        @test components[1].planet_shapiro

        @test isa(components[2], SolarWindDispersion)
        # @test components[3].model == 0

        @test isa(components[3], DispersionTaylor)

        @test isa(components[4], ExclusiveDispersionJump)

        @test isa(components[5], Spindown)

        @test isa(components[6], PhaseOffset)

        @test isa(components[7], DispersionMeasurementNoise)
    end

    @testset "form_residuals" begin
        params = model.param_handler._default_params_tuple
        wres = form_residuals(model, wtoas, params)
        @test all(abs(r[1]) < 4.1 * wtoa.toa.error for (r, wtoa) in zip(wres, wtoas))
        @test all(abs(r[2]) < 3.5 * wtoa.dminfo.error for (r, wtoa) in zip(wres, wtoas))
    end

    @testset "calc_chi2" begin
        calc_chi2_s = get_chi2_serial_func(model, wtoas)
        calc_chi2_p = get_chi2_parallel_func(model, wtoas)
        parv = read_param_values_to_vector(model.param_handler)
        # parnp = PyArray(parv)
        chi2_s = calc_chi2_s(parv)
        chi2_p = calc_chi2_p(parv)
        @test chi2_s ≈ chi2_p

        @test calc_chi2_reduced(model, wtoas, parv) < 1.5
    end

    @testset "calc_lnlike" begin
        calc_lnlike_s = get_lnlike_serial_func(model, wtoas)
        calc_lnlike_p = get_lnlike_parallel_func(model, wtoas)
        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model)
        # parnp = PyArray(parv)
        @test calc_lnlike_s(parv) ≈ calc_lnlike_p(parv)

        @test @ballocated(Vela.calc_lnlike_serial($model, $wtoas, $params)) == 0

        parv1 = read_param_values_to_vector(model.param_handler, params)
        parv1[end] *= 2
        @test calc_lnlike(model, wtoas, parv1) < calc_lnlike(model, wtoas, params)
    end

    @testset "priors" begin
        calc_lnprior = get_lnprior_func(model)
        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model.param_handler, params)
        @test isfinite(calc_lnprior(model.param_handler._default_params_tuple))
        @test calc_lnprior(params) == calc_lnprior(parv)

        calc_lnpost = get_lnpost_func(model, wtoas)
        @test isfinite(calc_lnpost(params))

        prior_transform = get_prior_transform_func(model)
        halfs = fill(0.5, length(parv))
        @test all(isfinite.(prior_transform(halfs)))
        # @test all(prior_transform(halfs) .≈ parv)
    end
end
