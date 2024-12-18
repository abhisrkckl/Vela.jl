# The `pyvela` command line interface

The `pyvela` script can be used to do simple analysis runs where fine control over
data handling, prior distributions, sampler, etc is not necessary. It has the following
syntax.

```
$ pyvela -h

A command line interface for the Vela.jl pulsar timing & noise analysis package

usage: pyvela [-h] [-P PRIOR_FILE] [--cheat_prior_scale CHEAT_PRIOR_SCALE] [-o OUTDIR] [-N NSTEPS] [-b BURNIN] [-t THIN] par_file tim_file

positional arguments:
  par_file
  tim_file

options:
  -h, --help            show this help message and exit
  -P PRIOR_FILE, --prior_file PRIOR_FILE
  --cheat_prior_scale CHEAT_PRIOR_SCALE
  -o OUTDIR, --outdir OUTDIR
  -N NSTEPS, --nsteps NSTEPS
  -b BURNIN, --burnin BURNIN
  -t THIN, --thin THIN
```
