"""Soil description classes."""
from __future__ import annotations

import dataclasses
from typing import Any, Callable, Generator, Optional, Sequence

import numpy as np
import scipy.stats as sst

from . import units as u
from ._typing import *

Q_ = u.Quantity

DispersivityFunction = Callable[["Soil", CharacteristicLength], Dispersivity]


def longitudinal_dispersivity_constant(
    dispersivity: Dispersivity,
) -> DispersivityFunction:
    """Longitudinal dispersivity of the soil.
    Based on a constant value.

    Args:
        dispersivity: Longitudinal dispersivity of the soil (dimensions [length]).

    Returns:
        Longitudinal dispersivity in same units as input L."""

    def _longitudinal_dispersivity_constant(
        soil: Soil, L: CharacteristicLength
    ) -> Dispersivity:
        """Longitudinal dispersivity of the soil.

        Args:
            soil (Soil): Soil for which the dispersivity is determined.
            L (u.QType[float]): Characteristic length scale of the soil (dimensions [length]).

        Returns:
            u.QType[float]: Longitudinal dispersivity in same units as input L."""
        return dispersivity.to(L.units)  # type: ignore

    return _longitudinal_dispersivity_constant


def longitudinal_dispersivity_Xu1995(
    soil: Soil, L: CharacteristicLength
) -> Dispersivity:
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
    return Dispersivity(Q_(0.83 * np.log10(L.m_as("m")) ** 2.414, "m").to(L.units))  # type: ignore


def relative_permeability(soil: Soil, theta: WaterContent) -> HydraulicConductivity:
    """Relative permeability of the soil (dimensionless)."""
    se = effective_saturation(soil, theta)
    l = soil.van_genuchten.l
    m = soil.van_genuchten.m

    return HydraulicConductivity(
        soil.K_sat * se**l * (1 - (1 - se ** (1 / m)) ** m) ** 2
    )


def _root_scalar(
    f: Callable[[float], float], bracket: Sequence[float], *args: Any, **kwargs: Any
) -> float:
    """Wrapper for scipy.optimize.root_scalar."""
    import scipy.optimize as opt  # type: ignore

    return opt.root_scalar(f, *args, bracket=bracket, **kwargs).root  # type: ignore


def volumetric_water_content(soil: Soil, q: Discharge) -> WaterContent:
    """Volumetric water content of the soil (dimensionless)."""
    import warnings

    if q > soil.K_sat:
        raise ValueError(
            "Groundwater recharge is too high for the soil under the unit gradient assumption"
        )

    # Under the assumption that dh/dz = 1, q = K_r
    K_r = q

    def rel_perm_error(theta: float) -> float:
        return (relative_permeability(soil, WaterContent(theta)) - K_r).m

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        try:
            # print(f"initial_theta = {initial_theta}")
            theta: float = _root_scalar(
                rel_perm_error,
                bracket=[soil.theta_r.m + 1e-9, soil.theta_s.m],
                method="brentq",
            )
        except RuntimeWarning as e:
            raise ValueError("Could not find a solution for the water content") from e

    return WaterContent(u.Q_(theta, "dimensionless"))


def effective_saturation(soil: Soil, theta: WaterContent) -> EffectiveSaturation:
    """Effective saturation of the soil (dimensionless)."""
    return EffectiveSaturation((theta - soil.theta_r) / (soil.theta_s - soil.theta_r))


def saturation(soil: Soil, theta: WaterContent) -> Saturation:
    """Saturation of the soil (dimensionless)."""
    return Saturation(theta / soil.porosity)


def bulk_density_Poelman1974(
    f_clay: FractionClay, f_oc: FractionOrganicCarbon, porosity: Porosity
) -> BulkDensity:
    """Bulk density estimate of the soil.

    Based on the method of Poelman (1974)

    Poelman, J. N. B. (1974). Particle density of river-clay soils
    http://edepot.wur.nl/110234

    Args:
        f_clay: Clay fraction of the soil (dimensionless).
        f_oc: Organic carbon fraction of the soil (dimensionless).
        theta_s: Saturated water content of the soil (dimensionless).

    Returns:
        Bulk density of the soil (dimensions [mass] / [volume])."""
    return BulkDensity(
        (1 - porosity).m_as("")
        * (f_oc * 1.47 + f_clay * 2.88 + (1 - f_oc - f_clay) * 2.66).m_as(""),
        "g/cm^3",
    )


def pressure_head(soil: Soil, theta: WaterContent) -> Q_[float]:
    """Pressure head for a given soil and volumetric water content

    Args:
        soil (pyfas.Soil): Soil
        theta (pint.Quantity[float]): Volumetric water content

    Returns:
        pint.Quantity[float]: Pressure head
    """
    S_e = effective_saturation(soil, theta)
    return (
        (S_e ** (1 / -soil.van_genuchten.m) - 1) ** (1 / soil.van_genuchten.n)  # type: ignore
    ) / soil.van_genuchten.alpha


def tortuosity_MillingtonQuirk1961(soil: Soil, theta: WaterContent) -> Tortuosity:
    """Tortuosity of the soil.
    Based on the empirical relationship of Millington and Quirk (1967).

    Millington, R. J., & Quirk, J. P. (1961). Permeability of porous solids.
        Transactions of the Faraday Society, 57, 1200. https://doi.org/10.1039/tf9615701200

    Args:
        soil: Soil for which the tortuosity is determined.
        theta: Saturated water content of the soil (dimensionless).

    Returns:
        Tortuosity of the soil (dimensionless)."""
    return Tortuosity((theta ** (7 / 3)) / (soil.theta_s**2))


@dataclasses.dataclass(frozen=True, eq=True)
class VanGenuchtenParameters:
    """Set of Van Genuchten Parameters"""

    alpha: Q_[float]
    """Van Genuchten parameter alpha (cm^-1)."""
    n: Q_[float]
    """Van Genuchten parameter n (dimensionless)."""
    l: Q_[float] = Q_(0.5, "dimensionless")
    """Van Genuchten parameter l (dimensionless)."""

    @property
    def m(self) -> Q_[float]:
        """Van Genuchten parameter m (dimensionless)."""
        return Q_(1 - 1 / self.n)


@dataclasses.dataclass(frozen=True, eq=True)
class TracerFitParameters:
    """Set of parameters for fitting the tracer-based Aaw."""

    x0: Q_[float]
    """x0 in cm^2/cm^3."""
    x1: Q_[float]
    """x1 in cm^2/cm^3."""
    x2: Q_[float]
    """x2 in cm^2/cm^3."""


@dataclasses.dataclass(frozen=True, eq=True)
class Soil:
    name: str = dataclasses.field(compare=True)
    """Name of the soil."""

    rho_b: BulkDensity = dataclasses.field(compare=False)
    """Bulk density of the soil in g/cm^3."""

    porosity: Porosity = dataclasses.field(compare=False)
    """Porosity of the soil (dimensionless)."""
    theta_s: SaturatedWaterContent = dataclasses.field(compare=False)
    """Saturated water content of the soil (dimensionless)."""
    theta_r: ResidualWaterContent = dataclasses.field(compare=False)
    """Residual water content of the soil (dimensionless)."""

    K_sat: SaturatedHydraulicConductivity = dataclasses.field(compare=False)
    """Saturated hydraulic conductivity of the soil in cm/s."""

    van_genuchten: VanGenuchtenParameters = dataclasses.field(compare=False)
    """Van Genuchten parameters."""

    f_oc: Optional[FractionOrganicCarbon] = dataclasses.field(
        compare=False, default=None
    )
    """Fraction of organic carbon in the soil (dimensionless)."""
    f_mo: Optional[FractionMetalOxides] = dataclasses.field(compare=False, default=None)
    """Fraction of metal oxides in the soil (dimensionless)."""
    f_sand: Optional[FractionSand] = dataclasses.field(compare=False, default=None)
    """Fraction of sand in the soil (dimensionless)."""
    f_clay: Optional[FractionClay] = dataclasses.field(compare=False, default=None)
    """Fraction of clay in the soil (dimensionless)."""
    f_silt: Optional[FractionSilt] = dataclasses.field(compare=False, default=None)
    """Fraction of silt in the soil (dimensionless)."""

    tracer_fit: Optional[TracerFitParameters] = dataclasses.field(
        compare=False, default=None
    )
    """Parameters for fitting the tracer-based Aaw."""
    soil_roughness_multiplier: Optional[Q_[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Scale factor for Aaw when using the thermodynamic method."""

    dispersivity: DispersivityFunction = dataclasses.field(
        default=longitudinal_dispersivity_Xu1995, compare=False
    )
    """Function for calculating the longitudinal dispersivity of the soil."""
    tortuosity: Callable[["Soil", WaterContent], Tortuosity] = dataclasses.field(
        default=tortuosity_MillingtonQuirk1961, compare=False
    )
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
        if self.soil_roughness_multiplier and self.soil_roughness_multiplier < 0:
            raise ValueError("Soil roughness multiplier must be positive.")
        if self.f_oc is not None and not (0 <= self.f_oc <= 1):
            raise ValueError("Fraction of organic carbon must be between 0 and 1.")
        if self.f_mo is not None and not (0 <= self.f_mo <= 1):
            raise ValueError("Fraction of metal oxides must be between 0 and 1.")
        if self.f_sand is not None and not (0 <= self.f_sand <= 1):
            raise ValueError("Fraction of sand must be between 0 and 1.")
        if self.f_clay is not None and not (0 <= self.f_clay <= 1):
            raise ValueError("Fraction of clay must be between 0 and 1.")
        if self.f_silt is not None and not (0 <= self.f_silt <= 1):
            raise ValueError("Fraction of silt must be between 0 and 1.")
        if (
            float(self.f_oc or 0.0)
            + float(self.f_mo or 0.0)
            + float(self.f_sand or 0.0)
            + float(self.f_clay or 0.0)
            + float(self.f_silt or 0.0)
            > 1
        ):
            raise ValueError(
                "Sum of fractions of organic carbon, metal oxides, sand, clay, and silt must be less than or equal to 1."
            )


def _theta_r(_: Porosity) -> ResidualWaterContent:
    return ResidualWaterContent(Q_(0.0))


def _theta_s(porosity: Porosity) -> SaturatedWaterContent:
    return SaturatedWaterContent(porosity)


@dataclasses.dataclass(frozen=True, eq=True)
class SoilClass:
    """Soil class with a name and a set of soils."""

    name: str = dataclasses.field(compare=True)
    """Name of the soil class."""

    f_clay: sst.rv_continuous = dataclasses.field(compare=False, repr=False)
    f_oc: sst.rv_continuous = dataclasses.field(compare=False, repr=False)
    f_sand: sst.rv_continuous = dataclasses.field(compare=False, repr=False)
    f_silt: sst.rv_continuous = dataclasses.field(compare=False, repr=False)
    f_mo: sst.rv_continuous = dataclasses.field(compare=False, repr=False)

    porosity: sst.rv_continuous = dataclasses.field(compare=False, repr=False)
    K_sat: sst.rv_continuous = dataclasses.field(compare=False, repr=False)

    van_genuchten_alpha: sst.rv_continuous = dataclasses.field(
        compare=False, repr=False
    )
    van_genuchten_n: sst.rv_continuous = dataclasses.field(compare=False, repr=False)
    van_genuchten_l: sst.rv_continuous = dataclasses.field(compare=False, repr=False)

    theta_r: Callable[[Porosity], ResidualWaterContent] = dataclasses.field(
        compare=False, repr=False, default=_theta_r
    )
    theta_s: Callable[[Porosity], SaturatedWaterContent] = dataclasses.field(
        compare=False, repr=False, default=_theta_s
    )

    rho_b: Callable[
        [FractionClay, FractionOrganicCarbon, Porosity], BulkDensity
    ] = dataclasses.field(compare=False, repr=False, default=bulk_density_Poelman1974)
    """Function for calculating the bulk density of the soil."""

    def mean(self) -> Soil:
        """Mean soil in the soil class."""
        porosity = Porosity(Q_(self.porosity.mean(), "dimensionless"))  # type: ignore
        f_c = FractionClay(Q_(self.f_clay.mean(), "dimensionless"))  # type: ignore
        f_oc = FractionOrganicCarbon(Q_(self.f_oc.mean(), "dimensionless"))  # type: ignore
        f_mo = FractionMetalOxides(Q_(self.f_mo.mean(), "dimensionless"))  # type: ignore
        f_sand = FractionSand(Q_(self.f_sand.mean(), "dimensionless"))  # type: ignore
        f_silt = FractionSilt(Q_(self.f_silt.mean(), "dimensionless"))  # type: ignore
        K_sat = SaturatedHydraulicConductivity(Q_(self.K_sat.mean(), "m/s"))  # type: ignore
        vg_alpha = Q_(float(self.van_genuchten_alpha.rvs(size=size, random_state=rng)), "1/cm")  # type: ignore
        vg_n = Q_(float(self.van_genuchten_n.rvs(size=size, random_state=rng)), "dimensionless")  # type: ignore
        vg_l = Q_(float(self.van_genuchten_l.rvs(size=size, random_state=rng)), "dimensionless")  # type: ignore
        return Soil(
            name=f"Mean of {self.__class__.__qualname__}(name={self.name})",
            rho_b=self.rho_b(f_c, f_oc, porosity),
            porosity=porosity,
            theta_s=self.theta_s(porosity),
            theta_r=self.theta_r(porosity),
            K_sat=K_sat,
            van_genuchten=VanGenuchtenParameters(alpha=vg_alpha, n=vg_n, l=vg_l),
            f_oc=f_oc,
            f_mo=f_mo,
            f_sand=f_sand,
            f_clay=f_c,
            f_silt=f_silt,
        )

    def sample(
        self, size: int = 1, rng: np.random.Generator = np.random.default_rng()
    ) -> Generator[Soil, None, None]:
        """Sample soils in the soil class."""
        _porosity = Porosity(Q_(self.porosity.rvs(size=size, random_state=rng), "dimensionless"))  # type: ignore
        _f_c = FractionClay(Q_(self.f_clay.rvs(size=size, random_state=rng), "dimensionless"))  # type: ignore
        _f_oc = FractionOrganicCarbon(Q_(self.f_oc.rvs(size=size, random_state=rng), "dimensionless"))  # type: ignore
        _f_mo = FractionMetalOxides(Q_(self.f_mo.rvs(size=size, random_state=rng), "dimensionless"))  # type: ignore
        _f_sand = FractionSand(Q_(self.f_sand.rvs(size=size, random_state=rng), "dimensionless"))  # type: ignore
        _f_silt = FractionSilt(Q_(self.f_silt.rvs(size=size, random_state=rng), "dimensionless"))  # type: ignore
        _K_sat = SaturatedHydraulicConductivity(Q_(self.K_sat.rvs(size=size, random_state=rng), "m/s"))  # type: ignore
        _vg_alpha = Q_(self.van_genuchten_alpha.rvs(size=size, random_state=rng), "1/cm")  # type: ignore
        _vg_n = Q_(self.van_genuchten_n.rvs(size=size, random_state=rng), "dimensionless")  # type: ignore
        _vg_l = Q_(self.van_genuchten_l.rvs(size=size, random_state=rng), "dimensionless")  # type: ignore

        for row in zip(
            _porosity,
            _f_c,
            _f_oc,
            _f_mo,
            _f_sand,
            _f_silt,
            _K_sat,
            _vg_alpha,
            _vg_n,
            _vg_l,
        ):
            porosity, f_c, f_oc, f_mo, f_sand, f_silt, K_sat, vg_alpha, vg_n, vg_l = row
            theta_r = self.theta_r(porosity)
            theta_s = self.theta_s(porosity)
            rho_b = self.rho_b(f_c, f_oc, porosity)
            yield Soil(
                name=f"Sample of {self.__class__.__qualname__}(name={self.name})",
                rho_b=rho_b,
                porosity=porosity,
                theta_s=theta_s,
                theta_r=theta_r,
                K_sat=K_sat,
                van_genuchten=VanGenuchtenParameters(alpha=vg_alpha, n=vg_n, l=vg_l),
                f_oc=f_oc,
                f_mo=f_mo,
                f_sand=f_sand,
                f_clay=f_c,
                f_silt=f_silt,
            )
