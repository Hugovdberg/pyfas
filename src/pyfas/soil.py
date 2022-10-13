import dataclasses
from typing import Callable, Optional

import numpy as np
import pint

from . import units as u

DispersivityFunction = Callable[["Soil", pint.Quantity[float]], pint.Quantity[float]]


def longitudinal_dispersivity_constant(
    dispersivity: pint.Quantity[float],
) -> DispersivityFunction:
    """Longitudinal dispersivity of the soil.
    Based on a constant value.

    Args:
        dispersivity: Longitudinal dispersivity of the soil (dimensions [length]).

    Returns:
        Longitudinal dispersivity in same units as input L."""

    def _longitudinal_dispersivity_constant(
        soil: "Soil", L: pint.Quantity[float]
    ) -> pint.Quantity[float]:
        """Longitudinal dispersivity of the soil.

        Args:
            soil (Soil): Soil for which the dispersivity is determined.
            L (pint.Quantity[float]): Characteristic length scale of the soil (dimensions [length]).

        Returns:
            pint.Quantity[float]: Longitudinal dispersivity in same units as input L."""
        return dispersivity.to(L.units)  # type: ignore

    return _longitudinal_dispersivity_constant


def longitudinal_dispersivity_Xu1995(
    soil: "Soil", L: pint.Quantity[float]
) -> pint.Quantity[float]:
    """Longitudinal dispersivity of the soil.
    Based on the empirical relationship of Xu and Eckstein (1995).

    Xu, M., Eckstein, Y., 1995. Use of weighted least-squares method in evaluation of the
        relationship between dispersivity and field scale. Groundwater 33, 905–908.
        https://doi.org/10.1111/j.1745-6584.1995.tb00035.x

    Args:
        soil: Soil for which the dispersivity is determined.
        L: Characteristic length scale of the soil (dimensions [length]).

    Returns:
        Longitudinal dispersivity in same units as input L."""
    return u.Q_(0.83 * np.log10(L.to("m").m) ** 2.414, "m").to(L.units)  # type: ignore


def relative_permeability(
    soil: "Soil", theta: pint.Quantity[float]
) -> pint.Quantity[float]:
    """Relative permeability of the soil (dimensionless)."""
    se = effective_saturation(soil, theta)
    l = soil.van_genuchten.l
    m = soil.van_genuchten.m

    return soil.K_sat * se**l * (1 - (1 - se ** (1 / m)) ** m) ** 2


def effective_saturation(
    soil: "Soil", theta: pint.Quantity[float]
) -> pint.Quantity[float]:
    """Effective saturation of the soil (dimensionless)."""
    return (theta - soil.theta_r) / (soil.theta_s - soil.theta_r)


def saturation(soil: "Soil", theta: pint.Quantity[float]) -> pint.Quantity[float]:
    """Saturation of the soil (dimensionless)."""
    return theta / soil.porosity


def tortuosity_MillingtonQuirk1961(
    soil: "Soil", theta: pint.Quantity[float]
) -> pint.Quantity[float]:
    """Tortuosity of the soil.
    Based on the empirical relationship of Millington and Quirk (1967).

    Millington, R. J., & Quirk, J. P. (1961). Permeability of porous solids.
        Transactions of the Faraday Society, 57, 1200. https://doi.org/10.1039/tf9615701200

    Args:
        soil: Soil for which the tortuosity is determined.
        theta: Saturated water content of the soil (dimensionless).

    Returns:
        Tortuosity of the soil (dimensionless)."""
    return (theta ** (7 / 3)) / (soil.theta_s**2)


@dataclasses.dataclass(frozen=True, eq=True)
class VanGenuchtenParameters:
    """Set of Van Genuchten Parameters"""

    alpha: pint.Quantity[float]
    """Van Genuchten parameter alpha (cm^-1)."""
    n: pint.Quantity[float]
    """Van Genuchten parameter n (dimensionless)."""
    l: pint.Quantity[float] = u.Q_(0.5, "dimensionless")
    """Van Genuchten parameter l (dimensionless)."""

    @property
    def m(self) -> pint.Quantity[float]:
        """Van Genuchten parameter m (dimensionless)."""
        return u.Q_(1 - 1 / self.n)


@dataclasses.dataclass(frozen=True, eq=True)
class TracerFitParameters:
    """Set of parameters for fitting the tracer-based Aaw."""

    x0: pint.Quantity[float]
    """x0 in cm^2/cm^3."""
    x1: pint.Quantity[float]
    """x1 in cm^2/cm^3."""
    x2: pint.Quantity[float]
    """x2 in cm^2/cm^3."""


@dataclasses.dataclass(frozen=True, eq=True)
class Soil:
    name: str = dataclasses.field(compare=True)
    """Name of the soil."""

    rho_b: pint.Quantity[float] = dataclasses.field(compare=False)
    """Bulk density of the soil in g/cm^3."""

    porosity: pint.Quantity[float] = dataclasses.field(compare=False)
    """Porosity of the soil (dimensionless)."""
    theta_s: pint.Quantity[float] = dataclasses.field(compare=False)
    """Saturated water content of the soil (dimensionless)."""
    theta_r: pint.Quantity[float] = dataclasses.field(compare=False)
    """Residual water content of the soil (dimensionless)."""

    K_sat: pint.Quantity[float] = dataclasses.field(compare=False)
    """Saturated hydraulic conductivity of the soil in cm/s."""

    van_genuchten: VanGenuchtenParameters = dataclasses.field(compare=False)
    """Van Genuchten parameters."""

    f_oc: Optional[pint.Quantity[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Fraction of organic carbon in the soil (dimensionless)."""
    f_mo: Optional[pint.Quantity[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Fraction of metal oxides in the soil (dimensionless)."""
    f_sand: Optional[pint.Quantity[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Fraction of sand in the soil (dimensionless)."""
    f_clay: Optional[pint.Quantity[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Fraction of clay in the soil (dimensionless)."""
    f_silt: Optional[pint.Quantity[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Fraction of silt in the soil (dimensionless)."""

    tracer_fit: Optional[TracerFitParameters] = dataclasses.field(
        compare=False, default=None
    )
    """Parameters for fitting the tracer-based Aaw."""
    soil_roughness_multiplier: Optional[pint.Quantity[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Scale factor for Aaw when using the thermodynamic method."""

    dispersivity: DispersivityFunction = dataclasses.field(
        default=longitudinal_dispersivity_Xu1995, compare=False
    )
    """Function for calculating the longitudinal dispersivity of the soil."""
    tortuosity: Callable[
        ["Soil", pint.Quantity[float]], pint.Quantity[float]
    ] = dataclasses.field(default=tortuosity_MillingtonQuirk1961, compare=False)
    """Function for calculating the tortuosity of the soil."""

    def __post_init__(self) -> None:
        if not (0 < self.porosity < 1):
            raise ValueError("Porosity must be between 0 and 1.")
        if not (0 <= self.theta_r <= self.porosity):
            raise ValueError(
                f"Residual water content must be between 0 and porosity ({self.porosity})."
            )
        if not (self.theta_r <= self.theta_s <= self.porosity):
            raise ValueError(
                f"Saturated water content must be between residual water content ({self.theta_r}) and porosity ({self.porosity})."
            )
