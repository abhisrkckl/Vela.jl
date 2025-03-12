@testset "WoodburyKernel" begin
    kernel = WoodburyKernel(
        WhiteNoiseKernel(),
        (
            PowerlawRedNoiseGP(10),
            PowerlawDispersionNoiseGP(15),
            PowerlawChromaticNoiseGP(5),
        ),
        randn(1000, 60),
    )
end
