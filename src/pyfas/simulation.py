import dataclasses
import datetime
import warnings
from typing import Any, Callable, Dict, Optional, ParamSpec, Tuple, TypeVar, cast

import numpy as np
import numpy.typing as npt

from . import _typing as t
from . import air_water_adsorption as awa_
from . import constants as c
from . import pfas as p
from . import soil as s
from . import solid_phase_adsorption as spa_
from . import units as u

Q_ = u.Quantity


def air_water_interfacial_area_thermodynamic(
    sim: "Simulation",
    theta: s.WaterContent,
) -> Q_[float]:
    """Calculate the air-water interfacial area."""
    return LTM(sim.soil, theta=theta, rho_w=sim.rho_w, sigma0=sim.sigma0)


def LTM(
    soil: s.Soil,
    *,
    theta: Q_[float],
    rho_w: Q_[float],
    sigma0: Q_[float],
) -> Q_[float]:
    if soil.soil_roughness_multiplier is None:
        raise ValueError(
            "soil_roughness_multiplier must be set for thermodynamic Aaw method."
        )

    g = Q_(9.81, "m/s^2")

    m = soil.van_genuchten.m
    alpha = soil.van_genuchten.alpha
    n = soil.van_genuchten.n
    Sr = (soil.theta_r / soil.theta_s).m_as("dimensionless")
    Sw = linspace((theta / soil.theta_s).m_as("dimensionless"), 1.0, 1000)

    def h(Sw: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return (((1 - Sr) / (Sw - Sr)) ** (1 / m) - 1) ** (1 / n) / alpha  # type: ignore

    Aaw = rho_w * g * soil.porosity / sigma0 * trapz(h(Sw), Sw)

    return Aaw * soil.soil_roughness_multiplier


def linspace(start: float, end: float, num: int) -> npt.NDArray[np.float64]:
    return np.linspace(start, end, num)  # type: ignore


def trapz(y: npt.NDArray[np.float64], x: npt.NDArray[np.float64]) -> float:
    return np.trapz(y, x)  # type: ignore


def air_water_interfacial_area_quadratic(
    sim: "Simulation",
    theta: s.WaterContent,
) -> Q_[float]:
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
    awa: awa_.AWAParameters

    name: str = dataclasses.field(compare=True)
    """Name of the simulation."""

    # T: Q_[float]
    # """Total simulation time in seconds."""
    L: s.CharacteristicLength
    """Characteristic length of the soil column for determining dispersion in cm."""

    ts: Q_[npt.NDArray[np.float64]]
    """Time steps for output in seconds."""
    zs: Q_[npt.NDArray[np.float64]]
    """Depth steps for output in cm."""
    Ci: Q_[npt.NDArray[np.float64]]
    """Initial concentration in mg/L."""

    q: t.Discharge
    """Infiltation rate in cm/s."""

    T0: Q_[float]
    """Duration of the contaminant pulse in s."""
    F0: Q_[float]
    """Total PFAS input during the contaminant pulse in µmol/cm^2."""

    C_rep: Q_[float] = Q_(0.0, "mg/L")
    """Representative aqueous concentration of PFAS in mg/L. If unknown, set to 0."""
    sigma0: Q_[float] = Q_(72.8, "mN/m")  # Water-air at 20 deg C
    """Surface tension of aqueous solution without PFAS in dyn/cm=mN/m."""
    rho_w: Q_[float] = Q_(1.0, "g/cm^3")
    """Density of water in g/cm^3."""
    Temp: Q_[float] = Q_(293.15, "K")
    """Temperature in K."""
    chi: c.Ionic = dataclasses.field(default=c.Ionic.NONE, compare=False)
    """Coefficient of ionization (chi) (dimensionless)."""

    air_water_interfacial_area: Callable[
        ["Simulation", s.WaterContent], Q_[float]
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

    theta: t.WaterContent
    """Water content in cm^3/cm^3."""

    Kaw: Q_[float]
    """Adsorption coefficient for air-water interfaces in cm^3/cm^2."""
    Aaw: Q_[float]
    """Air-water interfacial area in cm^2/cm^3."""
    Raw: Q_[float]
    """Retardation factor associated with air-water interfaces (dimensionless)."""

    Kd: Q_[float]
    """Distribution coefficient for adsorption to solids in cm^3/g."""
    Rs: Q_[float]
    """Retardation factor associated with adsorption to solids (dimensionless)."""

    R: Q_[float]
    """Total retardation factor (dimensionless)."""

    v: Q_[float]
    """Pore-flow velocity in cm/s."""

    D: Q_[float]
    """Dispersion coefficient in cm^2/s."""

    Ci: Q_[npt.NDArray[np.float64]]
    """Initial aqueous concentration in µmol/cm^3."""
    F0: Q_[float]
    """Total PFAS input during the contaminant pulse in µmol/cm^2."""

    raw_args: Tuple[Any, ...]
    """Raw arguments for the simulation."""
    raw_kwargs: Dict[str, Any]
    """Raw keyword arguments for the simulation."""
    raw_results: Optional[Tuple[Any, ...]]
    """Raw results from the simulation."""

    C_aq: Q_[npt.NDArray[np.float64]]
    """Aqueous concentration in µmol/cm^3 of pore water."""
    C_aq_bulk: Q_[npt.NDArray[np.float64]]
    """Aqueous concentration in µmol/cm^3 of porous medium."""
    C_aw: Q_[npt.NDArray[np.float64]]
    """Air-water concentration in µmol/cm^3 of porous medium."""
    C_s_eq: Q_[npt.NDArray[np.float64]]
    """Equilibrium solid-phase adsorption in µmol/cm^3 of porous medium."""
    C_s_kin: Q_[npt.NDArray[np.float64]]
    """Kinetic solid-phase adsorption in µmol/cm^3 of porous medium."""
    C_s: Q_[npt.NDArray[np.float64]]
    """Total solid-phase adsorption in µmol/cm^3 of porous medium."""
    C_total: Q_[npt.NDArray[np.float64]]
    """Total concentration in in µmol/cm^3 of porous medium."""

    start_time: datetime.datetime
    """Start time of the simulation."""
    end_time: datetime.datetime
    """End time of the simulation."""
    elapsed_time: datetime.timedelta
    """Elapsed time of the simulation."""


P = ParamSpec("P")
a = TypeVar("a")


def function_call(
    fun: Callable[P, a], *args: P.args, **kwargs: P.kwargs
) -> Tuple[Optional[a], Tuple[Tuple[Any, ...], Dict[str, Any]]]:
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

    theta = s.volumetric_water_content(sim.soil, sim.q)
    v = sim.q.to_velocity(theta)
    D = coefficient_of_dispersion(
        sim.pfas, sim.soil, sim.L, v, theta, sim.include_diffusion
    )
    return run_simulation(sim, theta, v, D)


def run_simulation(
    sim: Simulation, theta: t.WaterContent, v: t.Velocity, D: t.DispersionCoefficient
) -> SimulationResult:
    from pfas_leach_screening.analytical_soln import analytical_soln  # type: ignore

    Ci = u.try_convert_to_substance(sim.Ci, "µmol/cm^3", sim.pfas.M)
    C_rep = u.try_convert_to_substance(sim.C_rep, "µmol/cm^3", sim.pfas.M)

    match sim.Ci_method:
        case c.CiFlag.AQUEOUS:
            Kaw, Aaw, Raw, Kd, Rs, R_total = calculate_retardation(sim, theta, C_rep)
        case c.CiFlag.BULK:
            Kaw, Aaw, Raw, Kd, Rs, R_total = calculate_retardation(
                sim, theta, u.Q_(1e-10, "µmol/L")
            )
            if C_rep > 0:
                R_old = R_total
                for _ in range(10):
                    C_rep_ = C_rep / (R_total * theta)
                    Kaw, Aaw, Raw, Kd, Rs, R_total = calculate_retardation(
                        sim, theta, C_rep_
                    )
                    if (R_total - R_old).m < 1e-3:
                        break
                    R_old = R_total
                else:
                    raise RuntimeError(
                        "Failed to converge on C_rep within 10 iterations."
                    )

            Ci = Ci / (R_total * theta)

    pfas_total = u.try_convert_to_substance(sim.F0, "µmol/cm^2", sim.pfas.M)

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
        t=sim.ts.m_as("s"),
        z=sim.zs.m_as("cm"),
        t0=sim.T0.m_as("s"),
        pfas_tot=pfas_total.m_as("umol/cm^2"),
        Ci=Ci.m_as("umol/cm^3"),
        L=sim.L.m_as("cm"),
        v=v.m_as("cm/s"),
        theta=theta.m_as("dimensionless"),
        rhob=sim.soil.rho_b.m_as("g/cm^3"),
        D=D.m_as("cm^2/s"),
        Kd=Kd.m_as("cm^3/g"),
        alphas=alpha_s.m_as("1/s"),
        Fs=F_s.m_as("dimensionless"),
        R=R_total.m_as("dimensionless"),
        Raw=Raw.m_as("dimensionless"),
        Rs=Rs.m_as("dimensionless"),
        spaflag=sim.SPA_method,  # type: ignore
        cflag=sim.C_method,  # type: ignore
    )
    end_time = datetime.datetime.now()

    if result is None:
        (C_aq, C_s_mass, C_tot) = 3 * (
            cast(npt.NDArray[np.float64], np.array([np.nan], dtype=np.float64)),  # type: ignore
        )
    else:
        (C_aq, C_s_mass, C_tot) = result

    C_aq = Q_(C_aq, "umol/cm^3")
    C_s_mass = Q_(C_s_mass, "umol/g")
    C_tot = Q_(C_tot, "umol/cm^3")

    C_aq_bulk = C_aq * theta
    C_aw = C_aq * Kaw * Aaw
    C_s_eq = C_aq * F_s * Kd * sim.soil.rho_b
    C_s_kin = C_s_mass * sim.soil.rho_b
    C_s = C_s_eq + C_s_kin
    C_total = C_aq_bulk + C_aw + C_s_eq + C_s_kin

    if np.abs((C_total - C_tot).m_as("µmol/cm^3")).max() > 1e-6:  # type: ignore
        warnings.warn("Total concentration does not match sum of components.")

    return SimulationResult(
        simulation=sim,
        theta=theta,
        Kaw=Kaw,
        Aaw=Aaw,
        Raw=Raw,
        Kd=Kd,
        Rs=Rs,
        R=R_total,
        v=v,
        D=D,
        Ci=Ci,
        F0=pfas_total,
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


def calculate_retardation(
    sim: Simulation, theta: t.WaterContent, C_rep: Q_[float]
) -> Tuple[Q_[float], ...]:
    Kaw = sim.awa.Kaw(sim, C_rep)
    Aaw = sim.air_water_interfacial_area(sim, theta)
    Raw = Aaw * Kaw / theta

    Kd = sim.spa.Kd(C_rep)
    Rs = sim.soil.rho_b.q * Kd / theta.q

    R = 1 + Raw + Rs
    return Kaw, Aaw, Raw, Kd, Rs, R


def coefficient_of_dispersion(
    pfas: p.PFAS,
    soil: s.Soil,
    L: t.CharacteristicLength,
    v: t.Velocity,
    theta: s.WaterContent,
    include_diffusion: bool,
) -> t.DispersionCoefficient:
    if include_diffusion:
        if pfas.diffusivity is None:
            raise ValueError("Diffusion coefficient not provided.")
        diffusion = t.Diffusivity(soil.tortuosity(soil, theta) * pfas.diffusivity)
    else:
        diffusion = t.Diffusivity(0, "cm^2/s")
    return t.DispersionCoefficient(soil.dispersivity(soil, L).q * v + diffusion)
