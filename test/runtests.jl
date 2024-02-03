using Vela
using Test
using GeometricUnits
using LinearAlgebra
using Quadmath
using JuliaFormatter
using HDF5

const day_to_s = 86400

@testset "Vela" begin
    @testset "toa" begin
        ssb_obs_pos = distance.([18.0354099, 450.01472245, 195.05827732])
        ssb_obs_vel = speed.([-9.96231954e-05, 3.31555854e-06, 1.12968547e-06])
        obs_sun_pos = distance.([-15.89483533, -450.17185232, -195.18212616])
        obs_jupiter_pos = distance.([-1610.1796849, -1706.87348483, -681.22381513])
        obs_saturn_pos = distance.([-2392.85651431, 3109.13626083, 1405.71912274])
        obs_venus_pos = distance.([140.85922773, -217.65571843, -74.64804201])
        obs_uranus_pos = distance.([9936.62957939, -3089.07377113, -1486.17339104])
        obs_neptune_pos = distance.([11518.60924426, -9405.0693235, -4126.91030657])
        obs_earth_pos = distance.([0.01199435, 0.01159591, -0.01316261])

        @test_throws AssertionError EphemerisVectors(
            ssb_obs_vel,
            ssb_obs_vel,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )
        @test_throws AssertionError EphemerisVectors(
            ssb_obs_pos,
            ssb_obs_pos,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )
        @test_throws AssertionError EphemerisVectors(
            ssb_obs_pos,
            ssb_obs_vel,
            ssb_obs_vel,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )
        @test_throws AssertionError EphemerisVectors(
            ssb_obs_pos,
            1e6 * ssb_obs_vel,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )

        ephem_vecs = EphemerisVectors(
            ssb_obs_pos,
            ssb_obs_vel,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )

        # Aphelion distance
        @test sqrt(dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos)) < distance(509.0)
        @test sqrt(dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos)) < distance(509.0)

        # Angle
        @test acos(
            -dot(ephem_vecs.ssb_obs_pos, ephem_vecs.obs_sun_pos) / sqrt(
                dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos) *
                dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos),
            ),
        ) < 0.01

        # Speed of light
        @test dot(ephem_vecs.ssb_obs_vel, ephem_vecs.ssb_obs_vel) < 1

        toaval = time(parse(Float128, "4610197611.8464445127"))
        toaerr = time(1e-6)
        freq = frequency(1.4e9)
        phase = dimensionless(Float128(1000.0))

        @test_throws MethodError TOA(time(4610197611.8), toaerr, freq, phase)
        @test_throws AssertionError TOA(
            dimensionless(parse(Float128, "4610197611.8464445127")),
            toaerr,
            freq,
            phase,
        )
        @test_throws AssertionError TOA(toaval, dimensionless(1e-6), freq, phase)
        @test_throws AssertionError TOA(toaval, toaerr, time(1.4e9), phase)
        @test_throws MethodError TOA(toaval, toaerr, freq, dimensionless(1000.0))
        @test_throws AssertionError TOA(toaval, toaerr, freq, time(Float128(1000.0)))

        toa1 = TOA(toaval, toaerr, freq, phase)
        @test is_barycentered(toa1)
        @test !toa1.tzr
        @test toa1.level == 0

        toa2 = TOA(toaval, toaerr, freq, phase, ephem_vecs)
        @test !is_barycentered(toa2)
        @test !toa2.tzr
        @test toa2.level == 0

        dt = time(1.0)
        toa3 = correct_toa_delay(toa1, dt)
        @test toa3.value == toa1.value - dt
        @test toa3.error == toa1.error
        @test toa3.frequency == toa1.frequency
        @test toa3.phase == toa1.phase
        @test is_barycentered(toa3) == is_barycentered(toa1)
        @test toa3.tzr == toa1.tzr
        @test toa3.level == 1

        dphi = dimensionless(0.3)
        toa4 = correct_toa_phase(toa3, dphi)
        @test toa4.value == toa3.value
        @test toa4.error == toa3.error
        @test toa4.frequency == toa3.frequency
        @test toa4.phase == toa3.phase + dphi
        @test is_barycentered(toa4) == is_barycentered(toa3)
        @test toa4.tzr == toa3.tzr
        @test toa4.level == 2

        @testset "tzr_toa" begin
            tzrtoa = make_tzr_toa(toaval, freq, ephem_vecs)
            @test tzrtoa.tzr
            @test !is_barycentered(tzrtoa)
            @test tzrtoa.error == time(0.0)
            @test tzrtoa.level == 0
        end
    end

    @testset "read_model_and_toas" begin
        h5open("NGC6440E.hdf5") do f

            @testset "read_toas" begin
                toas = read_toas(f)
                @test !any([toa.tzr for toa in toas])
                @test length(toas) == 62
                @test all([toa.level == 0 for toa in toas])
                @test all([
                    frequency(1e9) < toa.frequency < frequency(2.5e9) for toa in toas
                ])
                @test all([
                    time(53470.0 * day_to_s) < toa.value < time(54200.0 * day_to_s) for
                    toa in toas
                ])
                @test all([modf(toa.phase.x)[1] == 0 for toa in toas])
                @test all([toa.error > time(0.0) for toa in toas])

            end

            @testset "read_tzr_toa" begin
                tzrtoa = read_tzr_toa(f)
                @test tzrtoa.tzr
                @test tzrtoa.level == 0
                @test tzrtoa.error > time(0.0)
                @test modf(tzrtoa.phase.x)[1] == 0
                @test frequency(1e9) < tzrtoa.frequency < frequency(2.5e9)
                @test time(53470.0 * day_to_s) < tzrtoa.value < time(54200.0 * day_to_s)
            end

            @testset "read_param_handler" begin
                param_handler = read_param_handler(f)
                @assert length(param_handler.params) >= length(param_handler._free_params)
                @assert length(param_handler.params) ==
                        length(param_handler._default_params_dict)
                @assert Set([p.name for p in param_handler._free_params]) == Set(["F0", "F1", "DECJ", "RAJ", "DM"])
            end

        end
    end

    # 
    #     @testset "read_toas" begin
    #         get_model_and_toas = pyimport("pint.models" => "get_model_and_toas")
    #     end
    #     model, toas = read_model_and_toas("NGC6440E.par", "NGC6440E.tim")
    #     @test !isempty(toas)
    #     @test all([toa.value.d == 1 for toa in toas])
    #     @test all([toa.error.d == 1 for toa in toas])
    #     @test all([toa.frequency.d == -1 for toa in toas])
    #     @test all([toa.phase.d == 0 for toa in toas])
    #     @test all([toa.level == 0 for toa in toas])
    #     @test all([!is_barycentered(toa) for toa in toas])
    #     @test all([!toa.tzr for toa in toas])

    #     @test "F0" in keys(model.param_handler._default_params_dict)
    #     @test "PHOFF" in keys(model.param_handler._default_params_dict)
    #     @test "NE_SW" in keys(model.param_handler._default_params_dict)
    #     @test !("DM1" in keys(model.param_handler._default_params_dict))
    #     @test !any([p.name == "PEPOCH" for p in model.param_handler._free_params])
    #     @test any([p.name == "PEPOCH" for p in model.param_handler.params])

    #     model, toas = read_model_and_toas("pure_rotator.par", "pure_rotator.tim")
    #     @test !isempty(toas)
    #     @test all([is_barycentered(toa) for toa in toas])

    #     @test "F0" in keys(model.param_handler._default_params_dict)
    #     @test "PHOFF" in keys(model.param_handler._default_params_dict)
    #     @test !any([p.name == "PEPOCH" for p in model.param_handler._free_params])
    #     @test any([p.name == "PEPOCH" for p in model.param_handler.params])
    # end

    @testset "parameter" begin
        @test_throws AssertionError Parameter(
            "PEPOCH",
            time(58000.0 * day_to_s),
            time(55000.0 * day_to_s),
            time(57000.0 * day_to_s),
            true,
        )

        pepoch = Parameter(
            "pepoch",
            time(56000.0 * day_to_s),
            time(55000.0 * day_to_s),
            time(57000.0 * day_to_s),
            true,
        )
        @assert pepoch.name == "PEPOCH"
        @test_throws AssertionError read_param(pepoch, 56000.0 * day_to_s)

        f0 = Parameter("F0", frequency(100.0), frequency(0.0), frequency(Inf), false)
        @test_throws AssertionError read_param(f0, -100.0)
        @test read_param(f0, 101.0) == frequency(101.0)
    end

    @testset "param handler" begin
        pepoch = Parameter(
            "pepoch",
            time(56000.0 * day_to_s),
            time(55000.0 * day_to_s),
            time(57000.0 * day_to_s),
            true,
        )
        f0 = Parameter("F0", frequency(100.0), frequency(0.0), frequency(Inf), false)
        f1 = Parameter("F1", GQ(-1e-14, -2), GQ(-Inf, -2), GQ(Inf, -2), false)

        @test_throws AssertionError ParamHandler([pepoch, f0, f0])

        param_handler = ParamHandler([pepoch, f0, f1])
        @test param_handler._free_params == [f0, f1]
        @test keys(param_handler._default_params_dict) == Set(["PEPOCH", "F0", "F1"])

        params_dict = read_params(param_handler, [100.01, -1.01e-14])
        @test keys(params_dict) == Set(["PEPOCH", "F0", "F1"])
    end

    # @testset "timing model" begin
    #     pepoch = Parameter(
    #         "pepoch",
    #         time(56000.0 * day_to_s),
    #         time(55000.0 * day_to_s),
    #         time(57000.0 * day_to_s),
    #         true,
    #     )
    #     f0 = Parameter("F0", frequency(100.0), frequency(0.0), frequency(Inf), false)
    #     f1 = Parameter("F1", GQ(-1e-14, -2), GQ(-Inf, -2), GQ(Inf, -2), false)

    #     param_handler = ParamHandler([pepoch, f0, f1])

    #     tzrtdb = time(Float128(55000.0))
    #     tzrfrq = frequency(1400.0)
    #     tzr_ephem = EphemerisVectors(
    #         distance.([18.0354099, 450.01472245, 195.05827732]),
    #         speed.([-9.96231954e-05, 3.31555854e-06, 1.12968547e-06]),
    #         distance.([-15.89483533, -450.17185232, -195.18212616]),
    #     )
    #     tzrtoa = make_tzr_toa(tzrtdb, tzrfrq, tzr_ephem)

    #     model = TimingModel(param_handler, tzrtoa)
    # end

    @testset "formatting" begin
        @test format(".")
    end
end
