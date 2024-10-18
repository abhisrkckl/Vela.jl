#!/usr/bin/env python

import sys
from pyvela.spnta import SPNTA

if __name__ == "__main__":
    parfile, timfile, jlsofile = sys.argv[1], sys.argv[2], sys.argv[3]
    SPNTA(parfile, timfile).save_jlso(jlsofile)
