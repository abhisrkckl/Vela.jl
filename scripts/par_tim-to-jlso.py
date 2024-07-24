#!/usr/bin/env python

from pint2vela import par_tim_to_jlso
import sys

parfile, timfile, jlsofile = sys.argv[1], sys.argv[2], sys.argv[3]

par_tim_to_jlso(parfile, timfile, jlsofile)
