"""Python interface for Vela.jl"""

from importlib import metadata

from .spnta import SPNTA
from .vela import vl as Vela

__version__ = metadata.version("pyvela")
