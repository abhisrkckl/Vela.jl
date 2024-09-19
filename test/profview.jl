using Vela
using ProfileView

const m, t = Vela.load_pulsar_data("datafiles/NGC6440E.jlso")
const lnlike = get_lnlike_serial_func(m, t)
const params = m.param_handler._default_params_tuple

function profile_test(n)
    for _ = 1:n
        _ = lnlike(params)
    end
end

# compilation
@profview profile_test(1)
# pure runtime
@profview profile_test(1000000)
