from typing import Dict

from . import constants as c
from . import pfas, soil
from . import solid_phase_adsorption as spa_
from . import registry as r

PFASs = {  # According to Guo et al. (2022), Table 2
    "PFPeA": pfas.PFAS(
        name="PFPeA",
        M=c.Q_(264.05, "g/mol"),
        Szyszkowski_a=c.Q_(3168.6, "mg/L") / c.Q_(264.05, "g/mol"),
        Szyszkowski_b=c.Q_(0.22),
        diffusion=c.Q_(12e-6, "cm**2/s"),
    ),
    "PFHxA": pfas.PFAS(
        name="PFHxA",
        M=c.Q_(314.05, "g/mol"),
        Szyszkowski_a=c.Q_(1350.42, "mg/L") / c.Q_(314.05, "g/mol"),
        Szyszkowski_b=c.Q_(0.21),
        diffusion=c.Q_(7.8e-6, "cm**2/s"),
    ),
    "PFHpA": pfas.PFAS(
        name="PFHpA",
        M=c.Q_(364.06, "g/mol"),
        Szyszkowski_a=c.Q_(345.86, "mg/L") / c.Q_(364.06, "g/mol"),
        Szyszkowski_b=c.Q_(0.22),
        diffusion=c.Q_(9.3e-6, "cm**2/s"),
    ),
    "PFOA": pfas.PFAS(
        name="PFOA",
        M=c.Q_(414.07, "g/mol"),
        Szyszkowski_a=c.Q_(62.11, "mg/L") / c.Q_(414.07, "g/mol"),
        Szyszkowski_b=c.Q_(0.19),
        diffusion=c.Q_(4.9e-6, "cm**2/s"),
    ),
    "PFNA": pfas.PFAS(
        name="PFNA",
        M=c.Q_(464.08, "g/mol"),
        Szyszkowski_a=c.Q_(5.11, "mg/L") / c.Q_(464.08, "g/mol"),
        Szyszkowski_b=c.Q_(0.16),
        diffusion=c.Q_(2.93e-6, "cm**2/s"),
    ),
    "PFDA": pfas.PFAS(
        name="PFDA",
        M=c.Q_(514.08, "g/mol"),
        Szyszkowski_a=c.Q_(3.7, "mg/L") / c.Q_(514.08, "g/mol"),
        Szyszkowski_b=c.Q_(0.17),
        diffusion=c.Q_(2.27e-6, "cm**2/s"),
    ),
    "PFBS": pfas.PFAS(
        name="PFBS",
        M=c.Q_(300.1, "g/mol"),
        Szyszkowski_a=c.Q_(2400.8, "mg/L") / c.Q_(300.1, "g/mol"),
        Szyszkowski_b=c.Q_(0.15),
        diffusion=c.Q_(11e-6, "cm**2/s"),
    ),
    "PFHxS": pfas.PFAS(
        name="PFHxS",
        M=c.Q_(400.12, "g/mol"),
        Szyszkowski_a=c.Q_(160.05, "mg/L") / c.Q_(400.12, "g/mol"),
        Szyszkowski_b=c.Q_(0.14),
        diffusion=c.Q_(4.5e-6, "cm**2/s"),
    ),
    "PFOS": pfas.PFAS(
        name="PFOS",
        M=c.Q_(500.13, "g/mol"),
        Szyszkowski_a=c.Q_(3.65, "mg/L") / c.Q_(500.13, "g/mol"),
        Szyszkowski_b=c.Q_(0.12),
        diffusion=c.Q_(5.4e-6, "cm**2/s"),
    ),
}
soils = {  # According to Guo et al. (2020), Table 1
    "Accusand": soil.Soil(
        name="Accusand",
        rho_b=c.Q_(1.65, "g/cm**3"),
        porosity=c.Q_(0.294),
        theta_s=c.Q_(0.294),
        theta_r=c.Q_(0.015),
        K_sat=c.Q_(2.0964e-2, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=c.Q_(0.04479, "1/cm"), n=c.Q_(4.0)
        ),
        tracer_fit=soil.TracerFitParameters(
            x0=c.Q_(633.96, "cm**2/cm**3"),
            x1=c.Q_(-1182.5, "cm**2/cm**3"),
            x2=c.Q_(548.54, "cm**2/cm**3"),
        ),
        soil_roughness_multiplier=c.Q_(4.15),
    ),
    "Vinton soil": soil.Soil(
        name="Vinton soil",
        rho_b=c.Q_(1.627, "g/cm**3"),
        porosity=c.Q_(0.395),
        theta_s=c.Q_(0.395),
        theta_r=c.Q_(0.056),
        K_sat=c.Q_(1.17e-3, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=c.Q_(0.02178, "1/cm"), n=c.Q_(3.451)
        ),
        tracer_fit=soil.TracerFitParameters(
            x0=c.Q_(1543.6, "cm**2/cm**3"),
            x1=c.Q_(-2848.6, "cm**2/cm**3"),
            x2=c.Q_(1305.0, "cm**2/cm**3"),
        ),
        soil_roughness_multiplier=c.Q_(4.15),
    ),
    "Schoonenburgse Heuvel - sand": soil.Soil(
        name="Schoonenburgse Heuvel - sand",
        rho_b=c.Q_(1.5, "g/cm**3"),
        porosity=c.Q_(0.427),
        theta_s=c.Q_(0.427),
        theta_r=c.Q_(0.02),
        K_sat=c.Q_(3.61e-4, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=c.Q_(0.0217, "1/cm"), n=c.Q_(1.735)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=c.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=c.Q_(5.0, "cm")
        ),
    ),
    "Schoonenburgse Heuvel - peat": soil.Soil(
        name="Schoonenburgse Heuvel - peat",
        rho_b=c.Q_(0.23, "g/cm**3"),
        porosity=c.Q_(0.85),
        theta_s=c.Q_(0.849),
        theta_r=c.Q_(0.01),
        K_sat=c.Q_(3.93519e-05, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=c.Q_(0.0119, "1/cm"), n=c.Q_(1.272)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=c.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=c.Q_(10.0, "cm")
        ),
    ),
    "Silva et al. (2020) - loam": soil.Soil(
        name="Silva et al. (2020) - loam",
        rho_b=c.Q_(1.33, "g/cm**3"),
        porosity=c.Q_(0.47),
        theta_s=c.Q_(0.43),
        theta_r=c.Q_(0.078),
        K_sat=c.Q_(2.89e-4, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=c.Q_(0.036, "1/cm"), n=c.Q_(1.56)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=c.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=c.Q_(35.0, "cm")
        ),
    ),
    "Silva et al. (2020) - Loamy sand": soil.Soil(
        name="Silva et al. (2020) - Loamy sand",
        rho_b=c.Q_(1.65, "g/cm**3"),
        porosity=c.Q_(0.44),
        theta_s=c.Q_(0.41),
        theta_r=c.Q_(0.057),
        K_sat=c.Q_(1.23e-3, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=c.Q_(0.0124, "1/cm"), n=c.Q_(2.28)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=c.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=c.Q_(35.0, "cm")
        ),
    ),
}

spa_matrix: Dict[
    soil.Soil, Dict[pfas.PFAS, spa_.SPAParameters]
] = {  # According to Guo et al. (2022), Table 1
    soils["Accusand"]: {
        PFASs["PFPeA"]: spa_.SPAParameters(
            pfas=PFASs["PFPeA"],
            soil=soils["Accusand"],
            Freundlich_K=c.Q_(0.0211, "(mg/kg)/(mg/L)**0.87"),
            Freundlich_N=c.Q_(0.87),
            frac_instant_adsorption=c.Q_(0.4),
            kinetic_adsorption_rate=c.Q_(5.9, "1/h"),
        ),
        PFASs["PFHxS"]: spa_.SPAParameters(
            pfas=PFASs["PFHxS"],
            soil=soils["Accusand"],
            Freundlich_K=c.Q_(0.0213, "(mg/kg)/(mg/L)**0.81"),
            Freundlich_N=c.Q_(0.81),
            frac_instant_adsorption=c.Q_(0.1),
            kinetic_adsorption_rate=c.Q_(3.1, "1/h"),
        ),
        PFASs["PFOA"]: spa_.SPAParameters(
            pfas=PFASs["PFOA"],
            soil=soils["Accusand"],
            Freundlich_K=c.Q_(0.1, "mg/kg/(mg/L)**0.87"),
            Freundlich_N=c.Q_(0.87),
            frac_instant_adsorption=c.Q_(0.4),
            kinetic_adsorption_rate=c.Q_(5.9, "1/h"),
        ),
        PFASs["PFOS"]: spa_.SPAParameters(
            pfas=PFASs["PFOS"],
            soil=soils["Accusand"],
            Freundlich_K=c.Q_(0.15, "mg/kg/(mg/L)**0.81"),
            Freundlich_N=c.Q_(0.81),
            frac_instant_adsorption=c.Q_(0.1),
            kinetic_adsorption_rate=c.Q_(3.1, "1/h"),
        ),
    },
    soils["Vinton soil"]: {
        PFASs["PFPeA"]: spa_.SPAParameters(
            pfas=PFASs["PFPeA"],
            soil=soils["Vinton soil"],
            Freundlich_K=c.Q_(0.122, "(mg/kg)/(mg/L)**0.87"),
            Freundlich_N=c.Q_(0.87),
            frac_instant_adsorption=c.Q_(0.16),
            kinetic_adsorption_rate=c.Q_(0.9, "1/h"),
        ),
        PFASs["PFHxS"]: spa_.SPAParameters(
            pfas=PFASs["PFHxS"],
            soil=soils["Vinton soil"],
            Freundlich_K=c.Q_(0.156, "(mg/kg)/(mg/L)**0.77"),
            Freundlich_N=c.Q_(0.77),
            frac_instant_adsorption=c.Q_(0.16),
            kinetic_adsorption_rate=c.Q_(0.9, "1/h"),
        ),
        PFASs["PFOA"]: spa_.SPAParameters(
            pfas=PFASs["PFOA"],
            soil=soils["Vinton soil"],
            Freundlich_K=c.Q_(0.58, "(mg/kg)/(mg/L)**0.87"),
            Freundlich_N=c.Q_(0.87),
            frac_instant_adsorption=c.Q_(0.16),
            kinetic_adsorption_rate=c.Q_(0.9, "1/h"),
        ),
        PFASs["PFOS"]: spa_.SPAParameters(
            pfas=PFASs["PFOS"],
            soil=soils["Vinton soil"],
            Freundlich_K=c.Q_(1.11, "(mg/kg)/(mg/L)**0.77"),
            Freundlich_N=c.Q_(0.77),
            frac_instant_adsorption=c.Q_(0.16),
            kinetic_adsorption_rate=c.Q_(0.9, "1/h"),
        ),
    },
}

PFASRegistry = r.PFASRegistry(
    pfas=PFASs,
    soil=soils,
    spa=spa_matrix,
)
