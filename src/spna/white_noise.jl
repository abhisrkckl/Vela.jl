function correct_resid(wn::MeasurementNoise, cres::ResidCorrectionBase, index)
    ef_idx = wn.efac_index_mask[index]
    efac = (ef_idx == 0) ? 1.0 : value(params.EFAC[ef_idx])

    eq_idx = wn.equad_index_mask[index]
    equad2 = (eq_idx == 0) ? 0.0 : value(params.EQUAD[eq_idx])

    return correct_resid(cres; efac = efac, equad2 = equad2)
end

function correct_resid(
    dmwn::DispersionMeasurementNoise,
    cres::WidebandResidCorrection,
    index,
)
    dmef_idx = dmwn.dmefac_index_mask[index]
    dmefac = (dmef_idx == 0) ? 1.0 : value(params.DMEFAC[dmef_idx])

    dmeq_idx = dmwn.dmequad_index_mask[index]
    dmequad2 = (dmeq_idx == 0) ? 0.0 : value(params.DMEQUAD[eq_idx])

    return correct_resid(cres; dmefac = dmefac, dmequad2 = dmequad2)
end
