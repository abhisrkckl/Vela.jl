export get_chi2_serial_func,
    get_chi2_parallel_func, get_lnlike_serial_func, get_lnlike_parallel_func

calc_chi2_serial(model::TimingModel, toas::Vector{TOA}, params::PyArray) =
    calc_chi2_serial(model, toas, Vector{Float64}(params))

get_chi2_serial_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_chi2_serial(model, toas, params)

calc_chi2(model::TimingModel, toas::Vector{TOA}, params::PyArray) =
    calc_chi2(model, toas, Vector{Float64}(params))

get_chi2_parallel_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_chi2(model, toas, params)

calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params::PyArray) =
    calc_lnlike_serial(model, toas, Vector{Float64}(params))

get_lnlike_serial_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_lnlike_serial(model, toas, params)

calc_lnlike(model::TimingModel, toas::Vector{TOA}, params::PyArray) =
    calc_lnlike(model, toas, Vector{Float64}(params))

get_lnlike_parallel_func(model, toas) = params -> calc_lnlike(model, toas, params)
