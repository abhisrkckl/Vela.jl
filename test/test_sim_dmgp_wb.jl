@testset "sim_dmgp_marg" begin
    model, toas = Vela.load_pulsar_data("datafiles/sim_dmgp_wb.jlso")

    pepoch_mjd = value(model.epoch) / day_to_s

    @testset "file save" begin
        Vela.save_pulsar_data("__test.jlso", model, toas)
        @test isfile("__test.jlso")
    end

    @testset "copy" begin
        m2 = TimingModel(
            model.pulsar_name,
            model.ephem,
            model.clock,
            model.units,
            model.epoch,
            model.components,
            model.kernel,
            model.param_handler,
            model.tzr_toa,
            model.priors,
        )
    end

    @testset "read_toas" begin
        @test !any([is_tzr(toa.toa) for toa in toas])
        @test length(toas) == 500
        @test all([
            frequency(0.49e9) < toa.toa.observing_frequency < frequency(1.51e9) for
            toa in toas
        ])
        @test all([
            time((50000.0 - pepoch_mjd) * day_to_s) <
            toa.toa.value <
            time((60000.1 - pepoch_mjd) * day_to_s) for toa in toas
        ])
        @test all([modf(toa.toa.pulse_number.x)[1] == 0 for toa in toas])
        @test all([toa.toa.error > time(0.0) for toa in toas])
    end

    @testset "param_handler" begin
        param_handler = model.param_handler
        @test Set(get_free_param_names(param_handler)) == Set([
            "RAJ",
            "DECJ",
            "PHOFF",
            "DM",
            "DM1",
            "F0",
            "F1",
            "TNDMAMP",
            "TNDMGAM",
            "TNREDAMP",
            "TNREDGAM",
            "DMEFAC1",
            "EFAC1",
        ])
        @test get_num_timing_params(model) == 7
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
        @test length(model.components) == 6
        @test !any(is_gp_noise, model.components)
    end

    @testset "kernel" begin
        @test isa(model.kernel, WoodburyKernel)
        @test isa(model.kernel.inner_kernel, WhiteNoiseKernel)
        @test length(model.kernel.gp_components) == 2
    end

    @testset "_calc_resids_and_Ndiag" begin
        params = model.param_handler._default_params_tuple
        y, Ndiag = Vela._calc_resids_and_Ndiag(model, toas, params)
        @test length(y) == length(Ndiag) == 2*length(toas)
        @test all(Ndiag .> 0)
    end

    @testset "calc_lnlike" begin
        calc_lnlike_s = get_lnlike_serial_func(model, toas)
        calc_lnlike_p = get_lnlike_parallel_func(model, toas)
        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model)
        # parnp = PyArray(parv)
        @test calc_lnlike_s(parv) ≈ calc_lnlike_p(parv)

        parv1 = read_param_values_to_vector(model.param_handler, params)
        parv1[end-2] *= 2
        @test calc_lnlike(model, toas, parv1) < calc_lnlike(model, toas, params)
    end
end
