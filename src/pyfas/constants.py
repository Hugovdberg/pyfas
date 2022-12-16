import enum


class AawFlag(enum.Flag):
    """Flag for the type of air-water interfacial area."""

    THERMODYNAMIC = 0
    TRACER = 1


class SPAFlag(enum.IntFlag):
    """Flag for the type of solid-phase adsorption."""

    EQUILIBRIUM = 0
    KINETIC = 1


class CFlag(enum.IntFlag):
    """Flag for the type of concentration calculation."""

    VOLUME_AVG = 0
    FLUX_AVG = 1


class CiFlag(enum.IntFlag):
    """Flag for the type of initial concentration."""

    AQUEOUS = 0
    BULK = 1


class Ionic(enum.IntFlag):
    """Flag for the type of ionic strength."""

    NONE = 1
    """Non-ionic PFAS or ionic PFAS in presence of swamping amount of electrolyte."""
    IONIC = 2
    """Ionic PFAS in presence of low electrolyte concentration."""
