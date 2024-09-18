@testset "components" begin
    include("test_solarsystem.jl")

    include("test_solarwind.jl")

    include("test_dispersion.jl")

    include("test_chromatic.jl")

    include("test_fd.jl")

    # include("test_wavex.jl")

    # include("test_ell1.jl")

    # include("test_dd.jl")

    # include("test_phoff.jl")

    # include("test_spindown.jl")

    # include("test_jump.jl")

    # include("test_glitch.jl")

    # include("test_white_noise.jl")
end
