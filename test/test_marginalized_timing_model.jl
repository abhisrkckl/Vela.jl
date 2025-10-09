@testset "MarginalizedTimingModel" begin
    weights = [1e40, 1e40, 1e40]
    pnames = ["PHOFF", "F0", "F1"]
    mtm = MarginalizedTimingModel(weights, pnames)
    @test is_gp_noise(mtm)
    @test length(get_marginalized_param_names(mtm)) == 3
    @test calc_noise_weights_inv(mtm, (;)) == (1 ./ weights)
    @test Vela.get_gp_npars(mtm) == 3
    print(mtm)
end
