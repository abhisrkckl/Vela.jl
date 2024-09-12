import numpy as np
from juliacall import Main as jl

jl.seval("import Vela")
jl.seval("using Distributions")
jl.seval("using DoubleFloats")
vl = jl.Vela


def to_jldd(x_: np.longdouble):
    """Convert a `np.longdouble` number to a `DoubleFloats.Double64` representation."""
    x = float(x_)
    y = float(x_ - x)
    return jl.Double64(x, y)
