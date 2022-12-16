"""Unit helper functions."""
from __future__ import annotations

import numbers
from typing import TYPE_CHECKING, Any, Generic, Self, Type, TypeVar, cast

import numpy as np
import numpy.typing as npt
import pint
import pint.facets.plain.definitions as pint_def

if TYPE_CHECKING:
    from pint.facets.plain.quantity import PlainQuantity as QType
    from pint.facets.plain.unit import PlainUnit as UType
else:
    from pint import Quantity as QType
    from pint import Unit as UType

a = TypeVar(
    "a", numbers.Number, float, int, npt.NDArray[np.float64], npt.NDArray[np.int64]
)
b = TypeVar(
    "b", numbers.Number, float, int, npt.NDArray[np.float64], npt.NDArray[np.int64]
)


def replace_percent(s: str) -> str:
    return str.replace(s, "%", " percent ")


ureg = pint.UnitRegistry(preprocessors=[replace_percent])
_percentage_unit = pint_def.UnitDefinition(
    name="percent",
    defined_symbol="%",
    aliases=tuple(),
    converter=pint_def.ScaleConverter(0.01),
    reference={},  # type: ignore
)
ureg.define(_percentage_unit)  # type: ignore


class Quantity(cast(Type[QType], ureg.Quantity), Generic[a]):  # type: ignore
    """Pint Quantity with some additional methods."""

    def __new__(
        cls, value: a | Quantity[a], units: UType | str | None = None
    ) -> Quantity[a]:
        return super().__new__(cls, value, units)  # type: ignore

    def to(
        self,
        other: UType | str | None = None,
        *contexts: str | pint.Context,
        **ctx_kwargs: Any,
    ) -> Self:
        """Return PlainQuantity rescaled to different units.

        Args:
            other: The unit to convert to. If None, the unit is not changed.
            contexts: Contexts to use for conversion.
            ctx_kwargs: Keyword arguments to pass to the contexts.

        Returns:
            The converted quantity.
        """
        return super().to(other, *contexts, **ctx_kwargs)  # type: ignore

    def m_as(
        self,
        units: UType | str,
    ) -> a:
        """Convert a quantity to a different unit and return the magnitude.

        Args:
            value: Quantity to convert.
            target_units: Target unit.

        Returns:
            Quantity in target unit.
        """
        return super().m_as(units)  # type: ignore

    @property
    def magnitude(self) -> a:
        return super().magnitude  # type: ignore

    @property
    def m(self) -> a:
        return super().m  # type: ignore

    def __mul__(self, other: Quantity[b] | b) -> Self:
        return super().__mul__(other)  # type: ignore

    __rmul__ = __mul__

    def __truediv__(self, other: Quantity[b] | b) -> Self:
        return super().__truediv__(other)  # type: ignore

    def __rtruediv__(self, other: Quantity[b] | b) -> Self:
        return super().__rtruediv__(other)  # type: ignore

    def __add__(self, other: Quantity[b] | b) -> Self:
        return super().__add__(other)  # type: ignore

    __radd__ = __add__

    def __sub__(self, other: Quantity[b] | b) -> Self:
        return super().__sub__(other)  # type: ignore

    def __rsub__(self, other: Quantity[b] | b) -> Self:
        return super().__rsub__(other)  # type: ignore

    def __pow__(self, other: Quantity[b] | b) -> Self:
        return super().__pow__(other)  # type: ignore


def try_convert_to_mass(
    value: Quantity[a],
    target_units: UType | str,
    molar_mass: Quantity[float],
) -> Quantity[a]:
    try:
        return value.to(target_units)
    except pint.DimensionalityError:
        mass_value = value * molar_mass
        return mass_value.to(target_units)


def try_convert_to_substance(
    value: Quantity[a],
    target_units: UType | str,
    molar_mass: Quantity[float],
) -> Quantity[a]:
    """Convert a concentration to µmol/cm^3.

    Args:
        C: Concentration in [substance]/[length]**3 or [mass]/[length]**3.
        M: Molecular weight of the substance in [mass]/[substance].

    Returns:
        Concentration in µmol/cm^3 if possible.
    """
    try:
        return value.to(target_units)
    except pint.DimensionalityError:
        substance_value = value / molar_mass
        return substance_value.to(target_units)


Q_ = Quantity
