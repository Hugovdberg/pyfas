"""API for accessing the various pyfas modules."""

from . import simulation, soil
from .constants import AawFlag, CFlag, CiFlag, Ionic, SPAFlag
from .data import PFASRegistry
from .pfas import PFAS
from .simulation import Simulation, SimulationResult, simulate
from .soil import Soil, TracerFitParameters, VanGenuchtenParameters
from .solid_phase_adsorption import SPAParameters
from .units import Q_, Quantity, units

__all__ = [
    "units",
    "AawFlag",
    "SPAFlag",
    "CFlag",
    "CiFlag",
    "Ionic",
    "Q_",
    "Quantity",
    "PFAS",
    "Soil",
    "TracerFitParameters",
    "VanGenuchtenParameters",
    "SPAParameters",
    "Simulation",
    "SimulationResult",
    "simulate",
    "PFASRegistry",
    "simulation",
    "soil",
]
