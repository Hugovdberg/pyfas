import abc
import dataclasses

import pint

from . import pfas as p
from . import soil as s
from . import units as u


@dataclasses.dataclass(frozen=True, eq=True)
class SPAParameters(abc.ABC):
    """Parameters for the solid-phase adsorption model."""

    pfas: p.PFAS = dataclasses.field(compare=True)
    """PFAS parameters."""

    soil: s.Soil = dataclasses.field(compare=True)
    """Soil parameters."""

    @abc.abstractmethod
    def Kd(self, C: pint.Quantity[float]) -> pint.Quantity[float]:
        """Equilibrium partition coefficient (Kd) in cm^3/g."""
        raise NotImplementedError


@dataclasses.dataclass(frozen=True, eq=True)
class KineticSorption:
    frac_instant_adsorption: pint.Quantity[float] = dataclasses.field(compare=False)
    """Fraction of PFAS that instantaneously adsorbs to the soil (dimensionless)."""

    kinetic_adsorption_rate: pint.Quantity[float] = dataclasses.field(compare=False)
    """Kinetic adsorption rate ($\\alpha_s$) in 1/s."""

    def __post_init__(self) -> None:
        if not (0 <= self.frac_instant_adsorption <= 1):
            raise ValueError(
                f"frac_instant_adsorption must be in the range [0, 1], but is {self.frac_instant_adsorption}."
            )


@dataclasses.dataclass(frozen=True, eq=True)
class FreundlichSorption(SPAParameters):

    Freundlich_K: pint.Quantity[float] = dataclasses.field(compare=False)
    """Freundlich isotherm constant in mg/kg/(mg/L)^N."""
    Freundlich_N: pint.Quantity[float] = dataclasses.field(compare=False)
    """Freundlich isotherm exponent (dimensionless)."""

    def Kd(self, C: pint.Quantity[float]) -> pint.Quantity[float]:
        return Freundlich_isotherm(self.Freundlich_K, self.Freundlich_N, C)


@dataclasses.dataclass(frozen=True, eq=True)
class KineticFreundlichSorption(FreundlichSorption, KineticSorption):
    pass


@dataclasses.dataclass(frozen=True, eq=True)
class FabregatPalauSorption(SPAParameters):
    def Kd(self, C: pint.Quantity[float]) -> pint.Quantity[float]:
        if self.soil.f_oc is None:
            raise ValueError(
                "The soil organic carbon fraction (f_oc) must be specified to use the Fabregat-Palau model."
            )
        if self.soil.f_silt is None:
            raise ValueError(
                "The soil silt fraction (f_silt) must be specified to use the Fabregat-Palau model."
            )
        if self.soil.f_clay is None:
            raise ValueError(
                "The soil clay fraction (f_clay) must be specified to use the Fabregat-Palau model."
            )
        K_oc = self.pfas.K_oc or K_oc_FabregatPalau2021(self.pfas.n_CF2)
        K_sc = self.pfas.K_sc or K_sc_FabregatPalau2021(self.pfas.n_CF2)

        return K_oc * self.soil.f_oc + K_sc * (self.soil.f_silt + self.soil.f_clay)  # type: ignore


@dataclasses.dataclass(frozen=True, eq=True)
class KineticFabregatPalauSorption(FabregatPalauSorption, KineticSorption):
    pass


@dataclasses.dataclass(frozen=True, eq=True)
class LinearSorption(SPAParameters):
    Kd_: pint.Quantity[float] = dataclasses.field(compare=True)
    """Equilibrium partition coefficient (Kd) in cm^3/g."""

    def Kd(self, C: pint.Quantity[float]) -> pint.Quantity[float]:
        return self.Kd_


@dataclasses.dataclass(frozen=True, eq=True)
class KineticLinearSorption(LinearSorption, KineticSorption):
    pass


def Freundlich_isotherm(
    Kf: pint.Quantity[float], N: pint.Quantity[float], C: pint.Quantity[float]
) -> pint.Quantity[float]:
    """Freundlich isotherm of the solid-phase adsorption coefficient.

    Args:
        Kf: Freundlich Kf [mg/kg / (mg/L)^N]
        N: Freundlich power
        C: Concentration in mg/L

    Returns:
        Solid-phase adsorption coefficient in cm^3/g."""
    return Kf * C ** (N - 1)  # type: ignore


def Kd_FabregatPalau(
    CF_2: pint.Quantity[int],
    f_oc: pint.Quantity[float],
    f_silt_clay: pint.Quantity[float],
) -> pint.Quantity[float]:
    K_oc = K_oc_FabregatPalau2021(CF_2)
    K_silt_clay = K_sc_FabregatPalau2021(CF_2)

    return K_oc * f_oc + K_silt_clay * f_silt_clay  # type: ignore


def K_sc_FabregatPalau2021(CF_2: pint.Quantity[int]) -> pint.Quantity[float]:
    return u.Q_(10 ** (0.32 * CF_2.m - 1.7), "L/kg")


def K_oc_FabregatPalau2021(CF_2: pint.Quantity[int]) -> pint.Quantity[float]:
    return u.Q_(10 ** (0.41 * CF_2.m - 0.7), "L/kg")
