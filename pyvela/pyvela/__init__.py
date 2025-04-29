"""Python interface for Vela.jl"""

from importlib import metadata

from .spnta import SPNTA
from .vela import vl as Vela

__version__ = metadata.version("pyvela")

assert (
    __version__ == Vela.pkg_version()
), f"Version mismatch between Vela and pyvela! Please check your installation. Vela version is {Vela.pkg_version()} and pyvela version is {__version__}."
