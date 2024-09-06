function save_pulsar_data(
    filename::String,
    model::TimingModel,
    toas::Vector{T},
) where {T<:TOABase}
    JLSO.save(filename, :model => model, :toas => toas)
end

function load_pulsar_data(filename::String)
    data = JLSO.load(filename)
    return data[:model], data[:toas]
end
