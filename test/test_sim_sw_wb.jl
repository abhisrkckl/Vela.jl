@testset "sim_sw.wb" begin
    model, wtoas = Vela.load_pulsar_data("datafiles/sim_sw.wb.jlso")

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
        @test !any([wtoa.toa.tzr for wtoa in wtoas])
        @test length(wtoas) == 500
        @test all([
            frequency(1.3e9) < wtoa.toa.observing_frequency < frequency(1.5e9) for
            wtoa in wtoas
        ])
        @test all([
            time(53999.0 * day_to_s) < wtoa.toa.value < time(56001.0 * day_to_s) for
            wtoa in wtoas
        ])
        @test all([modf(wtoa.toa.pulse_number.x)[1] == 0 for wtoa in wtoas])
        @test all([wtoa.toa.error > time(0.0) for wtoa in wtoas])
    end
end
