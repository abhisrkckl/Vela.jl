function save_pulsar_data(filename::String, model::TimingModel, toas::Vector{TOA})
    JLSO.save(filename, :model => model, :toas => toas)
end

function load_pulsar_data(filename::String)::Tuple{TimingModel,Vector{TOA}}
    data = JLSO.load(filename)
    return data[:model], data[:toas]
end
