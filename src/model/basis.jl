function basis_dot(basis, amplitudes, index)
    @assert length(amplitudes) == size(basis)[1]
    @assert index <= size(basis)[2]
    return dot(amplitudes, @view(basis[:, index]))
end
