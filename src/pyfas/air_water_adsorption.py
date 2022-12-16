import abc
import dataclasses
from typing import Protocol

import numpy as np

from . import constants as c
from . import pfas as p
from . import units as u

Q_ = u.Quantity


class _Simulation(Protocol):
    Temp: Q_[float]
    sigma0: Q_[float]
    chi: c.Ionic


@dataclasses.dataclass(frozen=True, eq=True)
class AWAParameters(abc.ABC):
    """Parameters for the solid-phase adsorption model."""

    @abc.abstractmethod
    def Kaw(self, sim: _Simulation, C: Q_[float]) -> Q_[float]:
        """Equilibrium partition coefficient (K_aw) in cm^3/cm^2."""
        raise NotImplementedError


@dataclasses.dataclass(frozen=True, eq=True)
class SzyszkowskiSorption(AWAParameters):
    """Szyszkowski sorption model parameters."""

    a: Q_[float] = dataclasses.field(compare=False)
    """Szyszkowski parameter a in µmol/cm^3."""
    b: Q_[float] = dataclasses.field(compare=False)
    """Szyszkowski parameter b (dimensionless)."""

    def Kaw(self, sim: _Simulation, C: Q_[float]) -> Q_[float]:
        """Equilibrium partition coefficient (K_aw) in cm^3/cm^2."""
        return Szyszkowski_isotherm(self, sim, C)


@dataclasses.dataclass(frozen=True, eq=True)
class LeSorption(AWAParameters):
    """Le sorption model parameters."""

    dG0: Q_[float] = dataclasses.field(compare=True)
    logKaw0: Q_[float] = dataclasses.field(compare=True)

    def Kaw(self, sim: _Simulation, C: Q_[float]) -> Q_[float]:
        """Equilibrium partition coefficient (K_aw) in cm^3/cm^2."""
        return Le_isotherm(self, sim, C)

    @classmethod
    def from_pfas(cls, pfas: p.PFAS) -> "LeSorption":
        """Create Le sorption parameters from PFAS parameters."""
        dG0 = Le_dG0(pfas)
        logKaw0 = Le_logKaw0(pfas)
        return cls(dG0, logKaw0)


def Szyszkowski_isotherm(
    awa: SzyszkowskiSorption, sim: _Simulation, C: Q_[float]
) -> Q_[float]:
    """Szyszkowski sorption model.

    Parameters
    ----------
    awa (AWAParameters): Air-Water Adsorption parameters.
    sim (_Simulation): Simulation parameters.
    C (u.QType[float]): Concentration of PFAS in water (µmol/L).

    Returns
    -------
    Kaw (u.QType[float]): Equilibrium partition coefficient (K_aw) in cm^3/cm^2.
    """
    R = Q_(1, "molar_gas_constant")
    K_aw = sim.sigma0 * awa.b / (int(sim.chi) * R * sim.Temp * (awa.a + C))
    return K_aw.to("cm^3/cm^2")


def Le_isotherm(awa: LeSorption, sim: _Simulation, C: Q_[float]) -> Q_[float]:
    """Le sorption model.

    Parameters
    ----------
    awa (AWAParameters): Air-Water Adsorption parameters.
    sim (_Simulation): Simulation parameters.
    C (u.QType[float]): Concentration of PFAS in water (µmol/L).

    Returns
    -------
    Kaw (u.QType[float]): Equilibrium partition coefficient (K_aw) in cm^3/cm^2.
    """
    R = Q_(1, "molar_gas_constant")
    omega = Q_(55.3, "mol/L")  # water molar density
    Keq: Q_[float] = 1 / omega * np.exp(-awa.dG0 / (R * sim.Temp.to("K")))  # type: ignore
    return Q_(10**awa.logKaw0.m, "cm^3/cm^2") / (1 + Keq * C)


def Le_dG0(pfas: p.PFAS) -> Q_[float]:
    """Le sorption model standard free energy of sorption.

    Parameters
    ----------
    pfas (p.PFAS): PFAS parameters.

    Returns
    -------
    dG0 (u.QType[float]): Standard free energy of sorption (kJ/mol).
    """
    kJ_per_mol = Q_(1.0, "kJ/mol")
    return (
        Q_(-14.29, "kJ/mol")
        + pfas.n_CFx * -3.57 * kJ_per_mol
        + pfas.n_CHx * -2.07 * kJ_per_mol
        + pfas.n_COO * 11.56 * kJ_per_mol
        + pfas.n_COOH * 0.34 * kJ_per_mol
        + pfas.n_SO3 * 11.48 * kJ_per_mol
        + pfas.n_R4N * 22.06 * kJ_per_mol
        + pfas.n_OH * 4.22 * kJ_per_mol
        + pfas.n_OSO3 * 10.78 * kJ_per_mol
        + pfas.n__O_ * 1.91 * kJ_per_mol
        + pfas.n__S_ * 1.79 * kJ_per_mol
        + pfas.n_NCH32CH2COO * 3.42 * kJ_per_mol
    )


def Le_logKaw0(pfas: p.PFAS) -> Q_[float]:
    """Le sorption model log of equilibrium partition coefficient at C = 0 mol/L.

    Parameters
    ----------
    pfas (p.PFAS): PFAS parameters.

    Returns
    -------
    logKaw0 (u.QType[float]): Log of equilibrium partition coefficient at C = 0 mol/L.
    """
    return (
        Q_(
            -5.19,
            "dimensionless",
        )
        + pfas.n_CFx * 0.60
        + pfas.n_CHx * 0.36
        + pfas.n_COO * -2.42
        + pfas.n_COOH * -0.47
        + pfas.n_SO3 * -2.35
        + pfas.n_R4N * -4.30
        + pfas.n_OH * -0.79
        + pfas.n_OSO3 * -2.39
        + pfas.n__O_ * -0.41
        + pfas.n__S_ * -0.21
        + pfas.n_NCH32CH2COO * -1.07
    )
