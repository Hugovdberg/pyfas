import dataclasses
from typing import Optional

import pint


@dataclasses.dataclass(frozen=True, eq=True)
class PFAS:
    name: str = dataclasses.field(compare=True)
    """Name of the PFAS."""

    M: pint.Quantity[float] = dataclasses.field(compare=False)
    """Molecular weight of PFAS in g/mol."""
    n_CF2: pint.Quantity[int] = dataclasses.field(compare=False)
    """Number of CF2 groups in PFAS (dimensionless)."""

    Szyszkowski_a: pint.Quantity[float] = dataclasses.field(compare=False)
    """Szyszkowski parameter a in µmol/cm^3."""
    Szyszkowski_b: pint.Quantity[float] = dataclasses.field(compare=False)
    """Szyszkowski parameter b (dimensionless)."""

    diffusion: pint.Quantity[float] = dataclasses.field(compare=False)
    """Diffusion coefficient (D0) in cm^2/s."""

    K_oc: Optional[pint.Quantity[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Adsorption coefficient on soil organic matter (Koc) in cm^3/g."""
    K_sc: Optional[pint.Quantity[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Adsorption coefficient on silt and clay (Ksc) in cm^3/g."""
    K_mo: Optional[pint.Quantity[float]] = dataclasses.field(
        compare=False, default=None
    )
    """Adsorption coefficient on metal oxides (Kmo) in cm^3/g."""
