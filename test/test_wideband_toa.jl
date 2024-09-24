@testset "wideband toa" begin
    toaval = time(parse(Double64, "4610197611.8464445127"))
    toaerr = time(1e-6)
    freq = frequency(1.4e9)
    pulse_number = dimensionless(Double64(1000.0))

    ephem = default_ephem()

    toa = TOA(toaval, toaerr, freq, pulse_number, ephem, 1)
    dminfo = DMInfo(GQ{-1}(1e16), GQ{-1}(1e11))
    wtoa = WidebandTOA(toa, dminfo)
    cwtoa1 = WidebandTOACorrection()

    @test dm_residual(wtoa.dminfo, cwtoa1.dm_correction) == dminfo.value
    @test scaled_dm_error_sqr(wtoa.dminfo, cwtoa1.dm_correction) == dminfo.error^Val(2)

    ddm = GQ{-1}(1e14)
    cdminfo = cwtoa1.dm_correction
    cdminfo2 = correct_dminfo(
        cdminfo;
        delta_dm = ddm,
        dmefac = dimensionless(1.1),
        dmequad2 = GQ{-2}(1e-12),
    )
    @test cdminfo2.model_dm == cdminfo.model_dm + ddm
    @test cdminfo2.dmefac == dimensionless(1.1)
    @test cdminfo2.dmequad2 == GQ{-2}(1e-12)
end
