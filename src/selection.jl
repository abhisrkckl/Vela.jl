export Selection, check

struct Selection
    mask::BitMatrix
end

check(sel::Selection, toa_index, param_index) = sel.mask[toa_index, param_index]
