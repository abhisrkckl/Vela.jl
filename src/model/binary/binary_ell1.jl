export BinaryELL1

"""A binary model representing a nearly circular orbit."""
struct BinaryELL1 <: BinaryELL1Base
    use_fbx::Bool
end

function ELL1State(ell1::BinaryELL1, ctoa::CorrectedTOA, params::NamedTuple)
    Δt = corrected_toa_value(ctoa) - params.TASC
    a1 = params.A1 + Δt * params.A1DOT
    ϵ1 = params.EPS1 + Δt * params.EPS1DOT
    ϵ2 = params.EPS2 + Δt * params.EPS2DOT

    Φ = mean_anomaly(Δt, params, ell1.use_fbx)
    Φ_trigs = sincos(Φ), sincos(2Φ), sincos(3Φ), sincos(4Φ)

    n = mean_motion(Δt, params, ell1.use_fbx)

    return ELL1State(Φ_trigs, n, a1, ϵ1, ϵ2)
end

"""Shapiro delay due to a nearly circular binary."""
function shapiro_delay(::BinaryELL1, state::ELL1State, params::NamedTuple)::GQ
    if !issubset((:M2, :SINI), keys(params))
        return time(0.0)
    end

    sinΦ = state.Φ_trigs[1][1]
    m2 = params.M2
    sini = params.SINI
    return -2 * m2 * log(1 - sini * sinΦ)
end

function show(io::IO, ell1::BinaryELL1)
    mode = ell1.use_fbx ? "FBX" : "PB"
    print(io, "BinaryELL1(with $mode)")
end
