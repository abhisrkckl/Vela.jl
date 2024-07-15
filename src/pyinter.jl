export get_chi2_serial_func,
    get_chi2_parallel_func, get_chi2_func, get_lnlike_serial_func, get_lnlike_parallel_func, get_lnlike_func

get_chi2_serial_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_chi2_serial(model, toas, params)

get_chi2_parallel_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_chi2(model, toas, params)

get_chi2_func(model, toas) =
    (nthreads() == 1) ? get_chi2_serial_func(model, toas) :
    get_chi2_parallel_func(model, toas)

get_lnlike_serial_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_lnlike_serial(model, toas, params)

get_lnlike_parallel_func(model, toas) = params -> calc_lnlike(model, toas, params)

get_lnlike_func(model, toas) =
    (nthreads() == 1) ? get_lnlike_serial_func(model, toas) :
    get_lnlike_parallel_func(model, toas)
