# Known issues

## Python multiprocessing and MPI

I have noticed segmentation faults while using `pyvela` with Python multiprocessing and/or MPI
even with Julia multithreading turned off. So parallelism doesn't work outside Julia.

