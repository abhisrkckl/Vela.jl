#!/usr/bin/env python

import sys

from pyvela.pint2vela import par_tim_to_jlso

parfile, timfile, jlsofile = sys.argv[1], sys.argv[2], sys.argv[3]

par_tim_to_jlso(parfile, timfile, jlsofile)
