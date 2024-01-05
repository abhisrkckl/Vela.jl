using GeometricUnits
using Quadmath
using PythonCall

export read_toas

function read_toas(parfile::String, timfile::String)
    get_model_and_toas = pyimport("pint.models" => "get_model_and_toas")
    setup_log = pyimport("pint.logging" => "setup")
    
    setup_log(level="WARNING")

    pint_model, pint_toas = get_model_and_toas(parfile, timfile)
    pint_toas.compute_pulse_numbers(pint_model)

    parse_quantity = (quantity_func, convert_func, values) -> [
        quantity_func(convert_func(val)) for val in values
    ]

    f128_convert = val -> parse(Float128, string(val))
    f64_convert = val -> pyconvert(Float64, val)
    vec3_convert = vec -> pyconvert(Vector{Float64}, vec)

    vals = parse_quantity(time, f128_convert, pint_toas.table["tdbld"].value * 86400)
    errs = parse_quantity(time, f64_convert, pint_toas.table["error"].quantity.si.value)
    freqs = parse_quantity(frequency, f64_convert, pint_toas.table["freq"].quantity.si.value)
    ssb_obs_poss = parse_quantity((xs) -> distance.(xs), vec3_convert, pint_toas.table["ssb_obs_pos"].quantity.to_value("lightsecond"))
    ssb_obs_vels = parse_quantity((xs) -> speed.(xs), vec3_convert, pint_toas.table["ssb_obs_vel"].quantity.to_value("lightsecond/second"))
    obs_sun_poss = parse_quantity((xs) -> distance.(xs), vec3_convert, pint_toas.table["obs_sun_pos"].quantity.to_value("lightsecond"))
    phases = parse_quantity(dimensionless, f128_convert, (pint_toas.table["delta_pulse_number"].value - pint_toas.table["pulse_number"].value))

    return [
        TOA(
            val, 
            err, 
            freq, 
            phase,
            iszero(ssb_obs_pos) ? nothing : EphemerisVectors(ssb_obs_pos, ssb_obs_vel, obs_sun_pos)
        )
        for (val, err, freq, phase, ssb_obs_pos, ssb_obs_vel, obs_sun_pos) in zip(vals, errs, freqs, phases, ssb_obs_poss, ssb_obs_vels, obs_sun_poss)
    ]
end