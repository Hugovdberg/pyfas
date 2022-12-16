import dataclasses
from typing import Optional

from . import units as u

Q_ = u.Quantity


@dataclasses.dataclass(frozen=True, eq=True)
class PFAS:
    name: str = dataclasses.field(compare=True)
    """Name of the PFAS."""

    M: Q_[float] = dataclasses.field(compare=False)
    """Molecular weight of PFAS in g/mol."""
    n_CFx: Q_[int] = dataclasses.field(compare=False)
    """Number of CF2 and CF3 groups in PFAS (dimensionless)."""

    # Szyszkowski_a: u.QType[float] = dataclasses.field(compare=False)
    # """Szyszkowski parameter a in µmol/cm^3."""
    # Szyszkowski_b: u.QType[float] = dataclasses.field(compare=False)
    # """Szyszkowski parameter b (dimensionless)."""

    diffusion: Q_[float] = dataclasses.field(compare=False)
    """Diffusion coefficient (D0) in cm^2/s."""

    K_oc: Optional[Q_[float]] = dataclasses.field(compare=False, default=None)
    """Adsorption coefficient on soil organic matter (Koc) in cm^3/g."""
    K_sc: Optional[Q_[float]] = dataclasses.field(compare=False, default=None)
    """Adsorption coefficient on silt and clay (Ksc) in cm^3/g."""
    K_mo: Optional[Q_[float]] = dataclasses.field(compare=False, default=None)
    """Adsorption coefficient on metal oxides (Kmo) in cm^3/g."""

    n_CHx: Q_[int] = dataclasses.field(compare=False, default=Q_(0, "dimensionless"))
    """Number of CH2 and CH3 groups in PFAS (dimensionless)."""
    n_COO: Q_[int] = dataclasses.field(compare=False, default=Q_(0, "dimensionless"))
    """Number of COO- groups in PFAS (dimensionless)."""
    n_COOH: Q_[int] = dataclasses.field(compare=False, default=Q_(0, "dimensionless"))
    """Number of COOH groups in PFAS (dimensionless)."""
    n_SO3: Q_[int] = dataclasses.field(compare=False, default=Q_(0, "dimensionless"))
    """Number of SO3- groups in PFAS (dimensionless)."""
    n_R4N: Q_[int] = dataclasses.field(compare=False, default=Q_(0, "dimensionless"))
    """Number of R4N+ groups in PFAS (dimensionless)."""
    n_OH: Q_[int] = dataclasses.field(compare=False, default=Q_(0, "dimensionless"))
    """Number of OH groups in PFAS (dimensionless)."""
    n_OSO3: Q_[int] = dataclasses.field(compare=False, default=Q_(0, "dimensionless"))
    """Number of OSO3- groups in PFAS (dimensionless)."""
    n__O_: Q_[int] = dataclasses.field(compare=False, default=Q_(0, "dimensionless"))
    """Number of -O- groups in PFAS (dimensionless)."""
    n__S_: Q_[int] = dataclasses.field(compare=False, default=Q_(0, "dimensionless"))
    """Number of -S- groups in PFAS (dimensionless)."""
    n_NCH32CH2COO: Q_[int] = dataclasses.field(
        compare=False, default=Q_(0, "dimensionless")
    )
    """Number of N(CH3)2-CH2-COO- groups in PFAS (dimensionless)."""


def estimate_diffusion_coefficient(n_CFx: int) -> Q_[float]:
    """Estimate the diffusion coefficient of a PFAS.

    Based on a linear regression of the diffusion coefficient of 10 PFAS against
    their structural properties.

    Parameters
    ----------
    n_CFx : int
        number of CF2 and CF3 groups in PFAS (dimensionless)

    Returns
    -------
    u.QType[float]
        Estimated diffusion coefficient in cm^2/s.
    """
    return Q_(10 ** (-4.5360 + -0.1088 * n_CFx), "cm^2/s")
