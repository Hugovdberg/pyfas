import dataclasses

import pint


@dataclasses.dataclass(frozen=True, eq=True)
class PFAS:
    name: str = dataclasses.field(compare=True)
    """Name of the PFAS."""

    M: pint.Quantity[float] = dataclasses.field(compare=False)
    """Molecular weight of PFAS in g/mol."""

    Szyszkowski_a: pint.Quantity[float] = dataclasses.field(compare=False)
    """Szyszkowski parameter a in µmol/cm^3."""
    Szyszkowski_b: pint.Quantity[float] = dataclasses.field(compare=False)
    """Szyszkowski parameter b (dimensionless)."""

    diffusion: pint.Quantity[float] = dataclasses.field(compare=False)
    """Diffusion coefficient (D0) in cm^2/s."""
