import dataclasses

import pint

from . import pfas as p
from . import soil as s


@dataclasses.dataclass(frozen=True, eq=True)
class SPAParameters:
    pfas: p.PFAS = dataclasses.field(compare=True)
    """PFAS parameters."""

    soil: s.Soil = dataclasses.field(compare=True)
    """Soil parameters."""

    Freundlich_K: pint.Quantity[float] = dataclasses.field(compare=False)
    """Freundlich isotherm constant in mg/kg/(mg/L)^N."""
    Freundlich_N: pint.Quantity[float] = dataclasses.field(compare=False)
    """Freundlich isotherm exponent (dimensionless)."""

    frac_instant_adsorption: pint.Quantity[float] = dataclasses.field(compare=False)
    """Fraction of PFAS that instantaneously adsorbs to the soil (dimensionless)."""

    kinetic_adsorption_rate: pint.Quantity[float] = dataclasses.field(compare=False)
    """Kinetic adsorption rate ($\alpha_s$) in 1/s."""

    def __post_init__(self) -> None:
        if not (0 <= self.frac_instant_adsorption <= 1):
            raise ValueError(
                f"frac_instant_adsorption must be in the range [0, 1], but is {self.frac_instant_adsorption}."
            )


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
    return Kf * C ** (N - 1)
