using PythonCall

function get_chi2_parallel_func(model, toas)
    return params -> begin
        result = -Inf

        # Release the GIL for parallel computation
        pythread = PythonCall.C.PyEval_SaveThread()
        try
            result = calc_chi2(model, toas, Vector{Float64}(params))
        finally
            PythonCall.C.PyEval_RestoreThread(pythread)
        end

        return result
    end
end

function get_lnlike_parallel_func(model, toas)
    return params -> begin
        result = -Inf

        # Release the GIL for parallel computation
        pythread = PythonCall.C.PyEval_SaveThread()
        try
            result = calc_lnlike(model, toas, Vector{Float64}(params))
        finally
            PythonCall.C.PyEval_RestoreThread(pythread)
        end

        return result
    end
end
