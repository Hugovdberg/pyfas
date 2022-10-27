"""API for accessing the various pyfas modules."""

from . import simulation, soil
from .constants import AawFlag, CFlag, CiFlag, Ionic, SPAFlag
from .data import PFASRegistry
from .pfas import PFAS
from .simulation import Simulation, SimulationResult, simulate
from .soil import Soil, TracerFitParameters, VanGenuchtenParameters
from .solid_phase_adsorption import (
    FabregatPalauSorption,
    FreundlichSorption,
    KineticFabregatPalauSorption,
    KineticFreundlichSorption,
    KineticLinearSorption,
    KineticSorption,
    LinearSorption,
    SPAParameters,
    K_oc_FabregatPalau2021,
    K_sc_FabregatPalau2021,
    Kd_FabregatPalau,
)
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
    "KineticSorption",
    "FabregatPalauSorption",
    "KineticFabregatPalauSorption",
    "FreundlichSorption",
    "KineticFreundlichSorption",
    "LinearSorption",
    "KineticLinearSorption",
    "K_oc_FabregatPalau2021",
    "K_sc_FabregatPalau2021",
    "Kd_FabregatPalau",
    "Simulation",
    "SimulationResult",
    "simulate",
    "PFASRegistry",
    "simulation",
    "soil",
]
