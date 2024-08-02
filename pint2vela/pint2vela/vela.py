from juliacall import Main as jl
import numpy as np

jl.seval("import Vela")
jl.seval("using Distributions")
jl.seval("using DoubleFloats")
vl = jl.Vela


def to_jldd(x_: np.longdouble):
    x = float(x_)
    y = float(x_ - x)
    return jl.Double64(x, y)
