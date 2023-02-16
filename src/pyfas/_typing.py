from . import units as u

Q_ = u.Quantity


class Porosity(u.Quantity[float]):
    """Porosity of a soil.

    Args:
        magnitude (float): magnitude of the porosity
        units (UnitLike, optional): units of the porosity, either as a pint.Unit or a string
    """

    _base_units = ""


class WaterContent(u.Quantity[float]):
    """Water content of a soil.

    Args:
        magnitude (float): magnitude of the water content
        units (UnitLike, optional): units of the water content, either as a pint.Unit or a string
    """

    _base_units = ""


class SaturatedWaterContent(u.Quantity[float]):
    """Saturated water content of a soil.

    Args:
        magnitude (float): magnitude of the saturated water content
        units (UnitLike, optional): units of the saturated water content, either as a pint.Unit or a string
    """

    _base_units = ""


class ResidualWaterContent(u.Quantity[float]):
    """Residual water content of a soil.

    Args:
        magnitude (float): magnitude of the residual water content
        units (UnitLike, optional): units of the residual water content, either as a pint.Unit or a string
    """

    _base_units = ""


class HydraulicConductivity(u.Quantity[float]):
    """Hydraulic conductivity of a soil.

    Args:
        magnitude (float): magnitude of the hydraulic conductivity
        units (UnitLike, optional): units of the hydraulic conductivity, either as a pint.Unit or a string
    """

    _base_units = "cm/day"


class SaturatedHydraulicConductivity(u.Quantity[float]):
    """Saturated hydraulic conductivity of a soil.

    Args:
        magnitude (float): magnitude of the saturated hydraulic conductivity
        units (UnitLike, optional): units of the saturated hydraulic conductivity, either as a pint.Unit or a string
    """

    _base_units = "cm/day"


class Saturation(u.Quantity[float]):
    """Saturation of a soil.

    Args:
        magnitude (float): magnitude of the saturation
        units (UnitLike, optional): units of the saturation, either as a pint.Unit or a string
    """

    _base_units = ""


class EffectiveSaturation(u.Quantity[float]):
    """Effective saturation of a soil.

    Args:
        magnitude (float): magnitude of the effective saturation
        units (UnitLike, optional): units of the effective saturation, either as a pint.Unit or a string
    """

    _base_units = ""


class CharacteristicLength(u.Quantity[float]):
    """Characteristic length of a soil.

    Args:
        magnitude (float): magnitude of the characteristic length
        units (UnitLike, optional): units of the characteristic length, either as a pint.Unit or a string
    """

    _base_units = "cm"


class Dispersivity(u.Quantity[float]):
    """Dispersivity of a soil.

    Args:
        magnitude (float): magnitude of the dispersivity
        units (UnitLike, optional): units of the dispersivity, either as a pint.Unit or a string
    """

    _base_units = "cm"


class DispersionCoefficient(u.Quantity[float]):
    """Dispersion coefficient of a soil.

    Args:
        magnitude (float): magnitude of the dispersion coefficient
        units (UnitLike, optional): units of the dispersion coefficient, either as a pint.Unit or a string
    """

    _base_units = "cm^2/day"


class Tortuosity(u.Quantity[float]):
    """Tortuosity of a soil.

    Args:
        magnitude (float): magnitude of the tortuosity
        units (UnitLike, optional): units of the tortuosity, either as a pint.Unit or a string
    """

    _base_units = ""


class FractionClay(u.Quantity[float]):
    """Fraction of clay in a soil.

    Args:
        magnitude (float): magnitude of the fraction of clay
        units (UnitLike, optional): units of the fraction of clay, either as a pint.Unit or a string
    """

    _base_units = ""


class FractionSilt(u.Quantity[float]):
    """Fraction of silt in a soil.

    Args:
        magnitude (float): magnitude of the fraction of silt
        units (UnitLike, optional): units of the fraction of silt, either as a pint.Unit or a string
    """

    _base_units = ""


class FractionSand(u.Quantity[float]):
    """Fraction of sand in a soil.

    Args:
        magnitude (float): magnitude of the fraction of sand
        units (UnitLike, optional): units of the fraction of sand, either as a pint.Unit or a string
    """

    _base_units = ""


class FractionOrganicCarbon(u.Quantity[float]):
    """Fraction of organic matter in a soil.

    Args:
        magnitude (float): magnitude of the fraction of organic matter
        units (UnitLike, optional): units of the fraction of organic matter, either as a pint.Unit or a string
    """

    _base_units = ""


class FractionMetalOxides(u.Quantity[float]):
    """Fraction of metal oxides in a soil.

    Args:
        magnitude (float): magnitude of the fraction of metal oxides
        units (UnitLike, optional): units of the fraction of metal oxides, either as a pint.Unit or a string
    """

    _base_units = ""


class BulkDensity(u.Quantity[float]):
    """Bulk density of a soil.

    Args:
        magnitude (float): magnitude of the bulk density
        units (UnitLike, optional): units of the bulk density, either as a pint.Unit or a string
    """

    _base_units = "kg/m^3"


class Velocity(u.Quantity[float]):
    """Average linear velocity through a soil.

    Args:
        magnitude (float): magnitude of the velocity
        units (UnitLike, optional): units of the velocity, either as a pint.Unit or a string
    """

    _base_units = "cm/day"

    def to_discharge(self, theta: WaterContent) -> "Discharge":
        """Convert the velocity to a discharge.

        Args:
            theta (WaterContent): water content of the soil

        Returns:
            Discharge: discharge through the soil
        """
        return Discharge(self * theta)


class Discharge(u.Quantity[float]):
    """Specific discharge through a soil.

    Args:
        magnitude (float): magnitude of the discharge
        units (UnitLike, optional): units of the discharge, either as a pint.Unit or a string
    """

    _base_units = "cm/day"

    def to_velocity(self, theta: WaterContent) -> Velocity:
        """Convert the discharge to a velocity.

        Args:
            theta (WaterContent): water content of the soil

        Returns:
            Velocity: average linear velocity through the soil
        """
        return Velocity(self / theta)


class Diffusivity(u.Quantity[float]):
    """Diffusivity of a pfas.

    Args:
        magnitude (float): magnitude of the diffusivity
        units (UnitLike, optional): units of the diffusivity, either as a pint.Unit or a string
    """

    _base_units = "cm^2/day"


__all__ = [
    "Porosity",
    "WaterContent",
    "SaturatedWaterContent",
    "ResidualWaterContent",
    "HydraulicConductivity",
    "SaturatedHydraulicConductivity",
    "Saturation",
    "EffectiveSaturation",
    "CharacteristicLength",
    "Dispersivity",
    "Tortuosity",
    "FractionClay",
    "FractionOrganicCarbon",
    "FractionSilt",
    "FractionSand",
    "FractionMetalOxides",
    "BulkDensity",
    "Velocity",
    "Discharge",
    "Diffusivity",
]
