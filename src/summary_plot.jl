using Vela
using GeometricUnits
using Plots

function plot_pulsar_summary(filename)
    model, toas = read_model_and_toas(filename)

    params = model.param_handler._default_params_tuple

    mjds = [Float64(value(toa.value)) for toa in toas] / day_to_s
    res = value.(Vela.form_residuals(model, toas, params))
    errs = [value(toa.error) for toa in toas]
    freqs = [value(toa.observing_frequency) for toa in toas]

    wres = res ./ errs

    p = scatter(
        mjds,
        res * 1e6,
        yerror = errs * 1e6,
        markershape = :cross,
        layout = (3, 1),
        plot_title = model.pulsar_name,
    )
    xlabel!(p, "MJD", subplot = 1)
    ylabel!(p, "Residuals (us)", subplot = 1)

    scatter!(
        p,
        freqs / 1e6,
        res * 1e6,
        yerror = errs * 1e6,
        markershape = :cross,
        subplot = 2,
    )
    xlabel!(p, "Frequency (MHz)", subplot = 2)
    ylabel!(p, "Residuals (us)", subplot = 2)

    histogram!(p, wres, subplot = 3)
    xlabel!(p, "Whitened residuals", subplot = 3)

    savefig("$(model.pulsar_name)_summary.pdf")
end
