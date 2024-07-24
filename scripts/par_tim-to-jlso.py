#!/usr/bin/env python

from pint2vela import read_model_and_toas, vl
from juliacall import Main as jl

import sys

jl.seval("import JLSO")

parfile, timfile, jlsofile = sys.argv[1], sys.argv[2], sys.argv[3]

mv, tv = read_model_and_toas(parfile, timfile)

jl.JLSO.save(
    jlsofile, 
    jl.Pair(
        jl.Symbol("model"),
        mv,
    ),
    jl.Pair(
        jl.Symbol("toas"),
        tv,
    )
)

