import dataclasses
import datetime
from typing import Any, Callable, Dict, Optional, ParamSpec, Tuple, TypeVar

import numpy as np
import numpy.typing as npt
import pint

from . import constants as c
from . import pfas as p
from . import soil as s
from . import solid_phase_adsorption as spa_
from . import units as u


def air_water_interfacial_area_thermodynamic(
    sim: "Simulation",
    theta: pint.Quantity[float],
) -> pint.Quantity[float]:
    """Calculate the air-water interfacial area."""
    soil = sim.soil
    if soil.soil_roughness_multiplier is None:
        raise ValueError(
            "soil_roughness_multiplier must be set for thermodynamic Aaw method."
        )

    g = u.Q_(9.81, "m/s^2")

    m = soil.van_genuchten.m
    alpha = soil.van_genuchten.alpha
    n = soil.van_genuchten.n
    Sr = soil.theta_r / soil.theta_s
    Sw = np.linspace(theta / soil.theta_s, 1.0, 1000)

    def h(Sw: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return (((1 - Sr) / (Sw - Sr)) ** (1 / m) - 1) ** (1 / n) / alpha  # type: ignore

    Aaw = sim.rho_w * g * soil.porosity / sim.sigma0 * np.trapz(h(Sw), Sw)

    return Aaw * soil.soil_roughness_multiplier


def air_water_interfacial_area_quadratic(
    sim: "Simulation",
    theta: pint.Quantity[float],
) -> pint.Quantity[float]:
    """Calculate the air-water interfacial area."""
    if sim.soil.tracer_fit is None:
        raise ValueError("Soil does not have a tracer fit.")
    Sw = s.saturation(soil=sim.soil, theta=theta)
    par = sim.soil.tracer_fit
    return par.x0 + par.x1 * Sw + par.x2 * Sw**2


@dataclasses.dataclass(frozen=True)
class Simulation:
    pfas: p.PFAS
    soil: s.Soil
    spa: spa_.SPAParameters
    name: str = dataclasses.field(compare=True)
    """Name of the simulation."""

    # T: pint.Quantity[float]
    # """Total simulation time in seconds."""
    L: pint.Quantity[float]
    """Length of the soil column in cm."""

    ts: pint.Quantity[npt.NDArray[np.float64]]
    """Time steps for output in seconds."""
    zs: pint.Quantity[npt.NDArray[np.float64]]
    """Depth steps for output in cm."""
    Ci: pint.Quantity[npt.NDArray[np.float64]]
    """Initial concentration in mg/L."""

    q: pint.Quantity[float]
    """Infiltation rate in cm/s."""

    T0: pint.Quantity[float]
    """Duration of the contaminant pulse in s."""
    F0: pint.Quantity[float]
    """Total PFAS input during the contaminant pulse in µmol/cm^2."""

    C_rep: pint.Quantity[float] = u.Q_(0.0, "mg/L")
    """Representative aqueous concentration of PFAS in mg/L. If unknown, set to 0."""
    sigma0: pint.Quantity[float] = u.Q_(72.8, "mN/m")  # Water-air at 20 deg C
    """Surface tension of aqueous solution without PFAS in dyn/cm=mN/m."""
    rho_w: pint.Quantity[float] = u.Q_(1.0, "g/cm^3")
    """Density of water in g/cm^3."""
    Temp: pint.Quantity[float] = u.Q_(293.15, "K")
    """Temperature in K."""
    chi: c.Ionic = dataclasses.field(default=c.Ionic.NONE, compare=False)
    """Coefficient of ionization (chi) (dimensionless)."""

    air_water_interfacial_area: Callable[
        ["Simulation", pint.Quantity[float]], pint.Quantity[float]
    ] = air_water_interfacial_area_thermodynamic
    """Flag for the type of air-water interfacial area."""
    Ci_method: c.CiFlag = c.CiFlag.AQUEOUS
    """Flag for the type of initial concentration."""
    SPA_method: c.SPAFlag = c.SPAFlag.KINETIC
    """Flag for the type of solid-phase adsorption."""
    C_method: c.CFlag = c.CFlag.VOLUME_AVG
    """Flag for the type of concentration calculation."""
    include_diffusion: bool = False
    """Flag for whether to include diffusion in the calculation of the dispersion."""


@dataclasses.dataclass(frozen=True)
class SimulationResult:
    """Simulation result with input and intermediary values"""

    simulation: Simulation
    """Simulation parameters."""

    theta: pint.Quantity[float]
    """Water content in cm^3/cm^3."""

    Kaw: pint.Quantity[float]
    """Adsorption coefficient for air-water interfaces in cm^3/cm^2."""
    Aaw: pint.Quantity[float]
    """Air-water interfacial area in cm^2/cm^3."""
    Raw: pint.Quantity[float]
    """Retardation factor associated with air-water interfaces (dimensionless)."""

    Kd: pint.Quantity[float]
    """Distribution coefficient for adsorption to solids in cm^3/g."""
    Rs: pint.Quantity[float]
    """Retardation factor associated with adsorption to solids (dimensionless)."""

    R: pint.Quantity[float]
    """Total retardation factor (dimensionless)."""

    v: pint.Quantity[float]
    """Pore-flow velocity in cm/s."""

    D: pint.Quantity[float]
    """Dispersion coefficient in cm^2/s."""

    Ci: pint.Quantity[npt.NDArray[np.float64]]
    """Initial aqueous concentration in µmol/cm^3."""
    F0: pint.Quantity[float]
    """Total PFAS input during the contaminant pulse in µmol/cm^2."""

    raw_args: Tuple[Any]
    """Raw arguments for the simulation."""
    raw_kwargs: Dict[str, Any]
    """Raw keyword arguments for the simulation."""
    raw_results: Optional[Tuple[Any, ...]]
    """Raw results from the simulation."""

    C_aq: pint.Quantity[npt.NDArray[np.float64]]
    """Aqueous concentration in µmol/cm^3 of pore water."""
    C_aq_bulk: pint.Quantity[npt.NDArray[np.float64]]
    """Aqueous concentration in µmol/cm^3 of porous medium."""
    C_aw: pint.Quantity[npt.NDArray[np.float64]]
    """Air-water concentration in µmol/cm^3 of porous medium."""
    C_s_eq: pint.Quantity[npt.NDArray[np.float64]]
    """Equilibrium solid-phase adsorption in µmol/cm^3 of porous medium."""
    C_s_kin: pint.Quantity[npt.NDArray[np.float64]]
    """Kinetic solid-phase adsorption in µmol/cm^3 of porous medium."""
    C_s: pint.Quantity[npt.NDArray[np.float64]]
    """Total solid-phase adsorption in µmol/cm^3 of porous medium."""
    C_total: pint.Quantity[npt.NDArray[np.float64]]
    """Total concentration in in µmol/cm^3 of porous medium."""

    start_time: datetime.datetime
    """Start time of the simulation."""
    end_time: datetime.datetime
    """End time of the simulation."""
    elapsed_time: datetime.timedelta
    """Elapsed time of the simulation."""


def volumetric_water_content(
    sim: Simulation,
) -> pint.Quantity[float]:
    """Volumetric water content of the soil (dimensionless)."""
    import scipy.optimize as opt
    import warnings

    # Under the assumption that dh/dz = 1, q = K_r
    K_r = sim.q

    def rel_perm_error(theta: float) -> float:
        return (s.relative_permeability(sim.soil, u.Q_(theta)) - K_r).m

    initial_theta = (sim.soil.theta_s.m + sim.soil.theta_r.m) / 2
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        theta: float = opt.fsolve(rel_perm_error, initial_theta)[0]
    return u.Q_(theta, "dimensionless")


def K_aw(sim: Simulation, C_rep: pint.Quantity[float]) -> pint.Quantity[float]:
    """Air-water interfacial adsorption coefficient [cm^3/cm^2].

    Args:
        sim (Simulation): Simulation parameters.

    Returns:
        float: Air-water interfacial adsorption coefficient [cm^3/cm^2]."""

    return (
        sim.sigma0
        * sim.pfas.Szyszkowski_b
        / (
            sim.chi
            * u.units.molar_gas_constant
            * sim.Temp
            * (sim.pfas.Szyszkowski_a + C_rep / sim.pfas.M)
        )
    )  # air-water interfacial coefficient [cm^3/cm^2]


P = ParamSpec("P")
a = TypeVar("a")


def function_call(
    fun: Callable[P, a], *args: P.args, **kwargs: P.kwargs
) -> Tuple[Optional[a], Tuple[Tuple[Any], Dict[str, Any]]]:
    """Call a function with the given arguments and return the result and the raw arguments."""
    try:
        return fun(*args, **kwargs), (args, kwargs)
    except Exception as e:
        print(f"Error: {e}")
        return None, (args, kwargs)


def simulate(
    sim: Simulation,
) -> SimulationResult:
    """Simulate PFAS transport in soil."""
    from pfas_leach_screening.analytical_soln import analytical_soln

    theta = volumetric_water_content(sim)
    v = (sim.q / theta).to("cm/s")
    D = coefficient_of_dispersion(sim, v, theta).to("cm^2/s")

    try:
        Ci = sim.Ci.to("umol/cm^3")
    except pint.errors.DimensionalityError:
        Ci = (sim.Ci / sim.pfas.M).to("umol/cm^3")

    match sim.Ci_method:
        case c.CiFlag.AQUEOUS:
            Kaw, Aaw, Raw, Kd, Rs, R = _calculate_retardation(sim, theta, sim.C_rep)
        case c.CiFlag.BULK:
            Kaw, Aaw, Raw, Kd, Rs, R = _calculate_retardation(
                sim, theta, u.Q_(1e-6, "mg/L")
            )
            if sim.C_rep > 0:
                R_old = R
                for _ in range(10):
                    C_rep = sim.C_rep / (R * theta)
                    Kaw, Aaw, Raw, Kd, Rs, R = _calculate_retardation(sim, theta, C_rep)
                    if (R - R_old).m < 1e-3:
                        break
                    R_old = R
                else:
                    raise RuntimeError(
                        "Failed to converge on C_rep within 10 iterations."
                    )

            Ci = Ci / (R * theta)

    try:
        F0 = sim.F0.to("umol/cm^2")
    except pint.errors.DimensionalityError:
        F0: pint.Quantity[float] = (sim.F0 / sim.pfas.M).to("umol/cm^2")

    match sim.SPA_method:
        case c.SPAFlag.KINETIC:
            if not isinstance(sim.spa, spa_.KineticSorption):
                raise TypeError(
                    "Kinetic sorption model selected but no kinetic sorption parameters provided."
                )
            F_s = sim.spa.frac_instant_adsorption
            alpha_s = sim.spa.kinetic_adsorption_rate
        case c.SPAFlag.EQUILIBRIUM:
            F_s = u.Q_(1.0, "dimensionless")
            alpha_s = u.Q_(0.0, "1/s")

    start_time = datetime.datetime.now()
    result, (args, kwargs) = function_call(
        analytical_soln,
        t=sim.ts.to("s").m,
        z=sim.zs.to("cm").m,
        t0=sim.T0.to("s").m,
        pfas_tot=F0.m,
        Ci=Ci.m,
        L=sim.L.to("cm").m,
        v=v.m,
        theta=theta.to("dimensionless").m,
        rhob=sim.soil.rho_b.to("g/cm^3").m,
        D=D.m,
        Kd=Kd.to("cm^3/g").m,
        alphas=alpha_s.to("1/s").m,
        Fs=F_s.to("dimensionless").m,
        R=R.to("dimensionless").m,
        Raw=Raw.to("dimensionless").m,
        Rs=Rs.to("dimensionless").m,
        spaflag=sim.SPA_method,  # type: ignore
        cflag=sim.C_method,  # type: ignore
    )
    end_time = datetime.datetime.now()

    if result is None:
        (C_aq, C_s_mass, C_tot) = 3 * (np.array([np.nan], dtype=np.float64),)
    else:
        (C_aq, C_s_mass, C_tot) = result

    C_aq = u.Q_(C_aq, "umol/cm^3")
    C_s_mass = u.Q_(C_s_mass, "umol/g")
    C_tot = u.Q_(C_tot, "umol/cm^3")

    C_aq_bulk = C_aq * theta
    C_aw = C_aq * Kaw * Aaw
    C_s_eq = C_aq * F_s * Kd * sim.soil.rho_b
    C_s_kin = C_s_mass * sim.soil.rho_b
    C_s = C_s_eq + C_s_kin
    C_total = C_aq_bulk + C_aw + C_s_eq + C_s_kin

    if (C_total - C_tot).max() > u.Q_(1e-6, "umol/cm^3"):
        print("Warning: C_total != C_tot")

    return SimulationResult(
        simulation=sim,
        theta=theta,
        Kaw=Kaw,
        Aaw=Aaw,
        Raw=Raw,
        Kd=Kd,
        Rs=Rs,
        R=R,
        v=v,
        D=D,
        Ci=Ci,
        F0=F0,
        C_aq=C_aq.to("umol/cm^3"),
        C_aq_bulk=C_aq_bulk.to("umol/cm^3"),
        C_aw=C_aw.to("umol/cm^3"),
        C_s_eq=C_s_eq.to("umol/cm^3"),
        C_s_kin=C_s_kin.to("umol/cm^3"),
        C_s=C_s.to("umol/cm^3"),
        C_total=C_total.to("umol/cm^3"),
        raw_args=args,
        raw_kwargs=kwargs,
        raw_results=result,
        start_time=start_time,
        end_time=end_time,
        elapsed_time=end_time - start_time,
    )


def _calculate_retardation(
    sim: Simulation, theta: pint.Quantity[float], C_rep: pint.Quantity[float]
) -> Tuple[pint.Quantity[float], ...]:
    Kaw = K_aw(sim, C_rep)
    Aaw = sim.air_water_interfacial_area(sim, theta)
    Raw = Aaw * Kaw / theta

    Kd = sim.spa.Kd(C_rep)
    Rs = sim.soil.rho_b * Kd / theta

    R = 1 + Raw + Rs
    return Kaw, Aaw, Raw, Kd, Rs, R


def coefficient_of_dispersion(
    sim: Simulation,
    v: pint.Quantity[float],
    theta: pint.Quantity[float],
) -> pint.Quantity[float]:
    if sim.include_diffusion:
        diffusion = sim.soil.tortuosity(sim.soil, theta) * sim.pfas.diffusion
    else:
        diffusion = 0
    return sim.soil.dispersivity(sim.soil, sim.L) * v + diffusion
