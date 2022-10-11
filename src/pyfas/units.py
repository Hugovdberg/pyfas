from typing import Generic

import pint
import pint._typing as pt
import pint.converters as pint_conv
import pint.definitions as pint_def

def replace_percent(s: str) -> str:
    return str.replace(s, "%", " percent ")

units = pint.UnitRegistry(preprocessors=[replace_percent])
units.define(
    pint_def.UnitDefinition("percent", "%", (), pint_conv.ScaleConverter(0.01))
)


class Quantity(units.Quantity, Generic[pt._MagnitudeType]):  # type: ignore
    pass


Q_ = Quantity
