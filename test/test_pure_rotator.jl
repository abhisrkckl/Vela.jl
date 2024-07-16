@testset "pure_rotator" begin
    model, toas = read_model_and_toas("pure_rotator.hdf5")

    @testset "model info" begin
        @test model.pulsar_name == "SIM0"
        @test model.ephem == "DE421"
        @test model.clock == "TT(BIPM2021)"
        @test model.units == "TDB"
    end

    @testset "read_toas" begin
        @test !any([toa.tzr for toa in toas])
        @test length(toas) == 100

        @test all([toa.observing_frequency == frequency(1.4e9) for toa in toas])
        @test all([
            time(53400.0 * day_to_s) < toa.value < time(56002.0 * day_to_s) for toa in toas
        ])
        @test all([modf(toa.pulse_number.x)[1] == 0 for toa in toas])
        @test all([toa.error > time(0.0) for toa in toas])
    end

    @testset "read_tzr_toa" begin
        tzrtoa = model.tzr_toa
        @test tzrtoa.tzr
        @test tzrtoa.error > time(0.0)
        @test modf(tzrtoa.pulse_number.x)[1] == 0
        @test frequency(1e9) < tzrtoa.observing_frequency < frequency(2.5e9)
        @test time(53400.0 * day_to_s) < tzrtoa.value < time(56002.0 * day_to_s)
    end

    @testset "read_param_handler" begin
        param_handler = model.param_handler
        @test Set(get_free_param_names(param_handler)) == Set(["F0", "F1", "PHOFF"])
        @test length(param_handler.multi_params) + length(param_handler.single_params) ==
              length(param_handler._default_params_tuple)
        @test length(get_free_param_names(param_handler)) ==
              length(param_handler._free_indices)
        @test sizeof(param_handler._default_params_tuple) ==
              sizeof(GQ{Float64}) * length(param_handler._default_quantities)

        params = model.param_handler._default_params_tuple
        parv = [params.PHOFF.x, params.F[1].x, params.F[2].x]
        @test read_param_values_to_vector(model.param_handler, params) == parv
    end

    @testset "read_components" begin
        components = model.components
        @test length(components) == 2
        @test isa(components[1], Spindown)
        @test isa(components[2], PhaseOffset)
    end

    @testset "repr" begin
        @test startswith(string(toas[1]), "TOA")
        display(toas)
        display(toas[1])
        @test startswith(string(model.tzr_toa), "TZRTOA")
        display(model.tzr_toa)
        @test startswith(string(model), "TimingModel")
        display(model)
    end

    params = model.param_handler._default_params_tuple
    parv = read_param_values_to_vector(model.param_handler, params)

    @testset "correct_toa" begin
        ctoa = correct_toa(model, toas[1], params)
        @test corrected_toa_value(ctoa) ≈ toas[1].value
        @test doppler_shifted_spin_frequency(ctoa) ≈ ctoa.spin_frequency
    end

    @testset "form_residual" begin
        resid = form_residual(model, toas[1], params, dimensionless(0.0))
        @test abs(resid) < time(1e-2)
    end

    @testset "calc_chi2" begin
        chi2 = calc_chi2(model, toas, params)
        nfree = length(get_free_param_names(model.param_handler))
        @test isa(chi2, Float64)
        @test chi2 / (length(toas) - nfree) < dimensionless(1.5)

        @test chi2 ≈ calc_chi2(model, toas, parv)
        @test chi2 ≈ Vela.calc_chi2_serial(model, toas, params)
        @test chi2 ≈ Vela.calc_chi2_serial(model, toas, parv)
        @test chi2 ≈ get_chi2_func(model, toas)(parv)

        @test @ballocated(Vela.calc_chi2_serial($model, $toas, $params)) == 0
    end

    @testset "calc_lnlike" begin
        lnlike_func = get_lnlike_parallel_func(model, toas)
        lnlike_serial_func = Vela.get_lnlike_serial_func(model, toas)

        # lnlike = calc_lnlike(model, toas, params)
        lnlike = lnlike_func(params)
        @test isa(lnlike, Float64)
        @test isfinite(lnlike)

        @test lnlike ≈ lnlike_func(parv)
        @test lnlike ≈ lnlike_serial_func(params)
        @test lnlike ≈ lnlike_serial_func(parv)
        @test lnlike ≈ get_lnlike_func(model, toas)(parv)

        @test @ballocated(Vela.calc_lnlike_serial($model, $toas, $params)) == 0
    end
end