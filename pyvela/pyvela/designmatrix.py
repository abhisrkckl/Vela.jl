from pint.models import TimingModel
from pint.toa import TOAs

def designmatrix(model: TimingModel, toas: TOAs):
    M_tm, M_tm_params, M_tm_units = model.designmatrix(toas)

    