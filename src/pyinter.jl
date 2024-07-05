export get_chi2_serial_func,
    get_chi2_parallel_func, get_lnlike_serial_func, get_lnlike_parallel_func

calc_chi2_serial(model::TimingModel, toas::Vector{TOA}, params::PyArray) =
    calc_chi2_serial(model, toas, Vector{Float64}(params))

get_chi2_serial_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_chi2_serial(model, toas, params)

calc_chi2(model::TimingModel, toas::Vector{TOA}, params::PyArray) =
    calc_chi2(model, toas, Vector{Float64}(params))

function get_chi2_parallel_func(model::TimingModel, toas::Vector{TOA})
    return params -> begin
        result = -Inf

        # Release the GIL for parallel computation
        pythread = PyEval_SaveThread()
        try
            result = calc_chi2(model, toas, params)
        finally
            PyEval_RestoreThread(pythread)
        end

        return result
    end
end

calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params::PyArray) =
    calc_lnlike_serial(model, toas, Vector{Float64}(params))

get_lnlike_serial_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_lnlike_serial(model, toas, params)

calc_lnlike(model::TimingModel, toas::Vector{TOA}, params::PyArray) =
    calc_lnlike(model, toas, Vector{Float64}(params))

function get_lnlike_parallel_func(model, toas)
    return params -> begin
        result = -Inf

        # Release the GIL for parallel computation
        pythread = PyEval_SaveThread()
        try
            result = calc_lnlike(model, toas, params)
        finally
            PyEval_RestoreThread(pythread)
        end

        return result
    end
end
