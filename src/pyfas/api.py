"""API for accessing the various pyfas modules."""

from . import simulation, soil
from .air_water_adsorption import (
    AWAParameters,
    Le_dG0,
    Le_isotherm,
    Le_logKaw0,
    LeSorption,
    Szyszkowski_isotherm,
    SzyszkowskiSorption,
)
from .constants import AawFlag, CFlag, CiFlag, Ionic, SPAFlag
from .data import PFASRegistry
from .pfas import PFAS
from .simulation import Simulation, SimulationResult, simulate
from .soil import Soil, TracerFitParameters, VanGenuchtenParameters
from .solid_phase_adsorption import (
    FabregatPalauSorption,
    FreundlichSorption,
    K_oc_FabregatPalau2021,
    K_sc_FabregatPalau2021,
    Kd_FabregatPalau,
    KineticFabregatPalauSorption,
    KineticFreundlichSorption,
    KineticLinearSorption,
    KineticSorption,
    LinearSorption,
    SPAParameters,
)
from .units import Q_, Quantity, ureg

__all__ = [
    "ureg",
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
    "AWAParameters",
    "SzyszkowskiSorption",
    "LeSorption",
    "Le_logKaw0",
    "Le_dG0",
    "Le_isotherm",
    "Szyszkowski_isotherm",
]
