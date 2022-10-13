from typing import Dict

from . import pfas
from . import registry as r
from . import soil
from . import solid_phase_adsorption as spa_
from . import units as u

PFASs = {
    "PFPeA": pfas.PFAS(
        name="PFPeA",
        M=u.Q_(264.05, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_a=u.Q_(3168.6, "mg/L")
        / u.Q_(264.05, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_b=u.Q_(0.22),  # According to Guo et al. (2022), Table 2
        diffusion=u.Q_(12e-6, "cm**2/s"),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(15, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(0.46, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CF2=u.Q_(4),
    ),
    "PFHxA": pfas.PFAS(
        name="PFHxA",
        M=u.Q_(314.05, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_a=u.Q_(1350.42, "mg/L")
        / u.Q_(314.05, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_b=u.Q_(0.21),  # According to Guo et al. (2022), Table 2
        diffusion=u.Q_(7.8e-6, "cm**2/s"),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(15, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(0.46, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CF2=u.Q_(5),
    ),
    "PFHpA": pfas.PFAS(
        name="PFHpA",
        M=u.Q_(364.06, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_a=u.Q_(345.86, "mg/L")
        / u.Q_(364.06, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_b=u.Q_(0.22),  # According to Guo et al. (2022), Table 2
        diffusion=u.Q_(9.3e-6, "cm**2/s"),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(50, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=None,
        K_mo=None,
        n_CF2=u.Q_(6),
    ),
    "PFOA": pfas.PFAS(
        name="PFOA",
        M=u.Q_(414.07, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_a=u.Q_(62.11, "mg/L")
        / u.Q_(414.07, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_b=u.Q_(0.19),  # According to Guo et al. (2022), Table 2
        diffusion=u.Q_(4.9e-6, "cm**2/s"),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(107, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(3.3, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CF2=u.Q_(7),
    ),
    "PFNA": pfas.PFAS(
        name="PFNA",
        M=u.Q_(464.08, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_a=u.Q_(5.11, "mg/L")
        / u.Q_(464.08, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_b=u.Q_(0.16),  # According to Guo et al. (2022), Table 2
        diffusion=u.Q_(2.93e-6, "cm**2/s"),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(324, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(2.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CF2=u.Q_(8),
    ),
    "PFDA": pfas.PFAS(
        name="PFDA",
        M=u.Q_(514.08, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_a=u.Q_(3.7, "mg/L")
        / u.Q_(514.08, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_b=u.Q_(0.17),  # According to Guo et al. (2022), Table 2
        diffusion=u.Q_(2.27e-6, "cm**2/s"),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(604.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(14.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CF2=u.Q_(9),
    ),
    "PFBS": pfas.PFAS(
        name="PFBS",
        M=u.Q_(300.1, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_a=u.Q_(2400.8, "mg/L")
        / u.Q_(300.1, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_b=u.Q_(0.15),  # According to Guo et al. (2022), Table 2
        diffusion=u.Q_(11e-6, "cm**2/s"),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(11.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(0.44, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CF2=u.Q_(4),
    ),
    "PFHxS": pfas.PFAS(
        name="PFHxS",
        M=u.Q_(400.12, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_a=u.Q_(160.05, "mg/L")
        / u.Q_(400.12, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_b=u.Q_(0.14),  # According to Guo et al. (2022), Table 2
        diffusion=u.Q_(4.5e-6, "cm**2/s"),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(50, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(1.2, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CF2=u.Q_(6),
    ),
    "PFOS": pfas.PFAS(
        name="PFOS",
        M=u.Q_(500.13, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_a=u.Q_(3.65, "mg/L")
        / u.Q_(500.13, "g/mol"),  # According to Guo et al. (2022), Table 2
        Szyszkowski_b=u.Q_(0.12),  # According to Guo et al. (2022), Table 2
        diffusion=u.Q_(5.4e-6, "cm**2/s"),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(609, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(9.4, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CF2=u.Q_(8),
    ),
}
soils = {  # According to Guo et al. (2020), Table 1
    "Accusand": soil.Soil(
        name="Accusand",
        rho_b=u.Q_(1.65, "g/cm**3"),
        porosity=u.Q_(0.294),
        theta_s=u.Q_(0.294),
        theta_r=u.Q_(0.015),
        K_sat=u.Q_(2.0964e-2, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.04479, "1/cm"), n=u.Q_(4.0)
        ),
        tracer_fit=soil.TracerFitParameters(
            x0=u.Q_(633.96, "cm**2/cm**3"),
            x1=u.Q_(-1182.5, "cm**2/cm**3"),
            x2=u.Q_(548.54, "cm**2/cm**3"),
        ),
        soil_roughness_multiplier=u.Q_(4.15),
        f_oc=u.Q_(0.04, "percent"),  # According to Guo et al. (2020), Section 4
        f_mo=u.Q_(14.0, "ug/g")  # FeOx
        + u.Q_(2.5, "ug/g")  # MnOx
        + u.Q_(12.0, "ug/g"),  # AlOx; According to Guo et al. (2020), Section 4
        f_clay=u.Q_(0.0, "percent"),  # According to Guo et al. (2020), Section 4
        f_silt=u.Q_(0.0, "percent"),  # According to Guo et al. (2020), Section 4
    ),
    "Vinton soil": soil.Soil(
        name="Vinton soil",
        rho_b=u.Q_(1.627, "g/cm**3"),
        porosity=u.Q_(0.395),
        theta_s=u.Q_(0.395),
        theta_r=u.Q_(0.056),
        K_sat=u.Q_(1.17e-3, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.02178, "1/cm"), n=u.Q_(3.451)
        ),
        tracer_fit=soil.TracerFitParameters(
            x0=u.Q_(1543.6, "cm**2/cm**3"),
            x1=u.Q_(-2848.6, "cm**2/cm**3"),
            x2=u.Q_(1305.0, "cm**2/cm**3"),
        ),
        soil_roughness_multiplier=u.Q_(4.15),
        f_oc=u.Q_(0.1, "percent"),  # According to Guo et al. (2020), Section 4
        f_mo=u.Q_(0.0, "percent"),  # According to Guo et al. (2020), Section 4
        f_clay=u.Q_(4.7, "percent"),  # According to Guo et al. (2020), Section 4
        f_silt=u.Q_(0.0, "percent"),  # According to Guo et al. (2020), Section 4
    ),
    "Schoonenburgse Heuvel - sand": soil.Soil(
        name="Schoonenburgse Heuvel - sand",
        rho_b=u.Q_(1.5, "g/cm**3"),
        porosity=u.Q_(0.427),
        theta_s=u.Q_(0.427),
        theta_r=u.Q_(0.02),
        K_sat=u.Q_(3.61e-4, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0217, "1/cm"), n=u.Q_(1.735)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=u.Q_(5.0, "cm")
        ),
    ),
    "Schoonenburgse Heuvel - peat": soil.Soil(
        name="Schoonenburgse Heuvel - peat",
        rho_b=u.Q_(0.23, "g/cm**3"),
        porosity=u.Q_(0.85),
        theta_s=u.Q_(0.849),
        theta_r=u.Q_(0.01),
        K_sat=u.Q_(3.93519e-05, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0119, "1/cm"), n=u.Q_(1.272)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=u.Q_(10.0, "cm")
        ),
    ),
    "Silva et al. (2020) - loam": soil.Soil(
        name="Silva et al. (2020) - loam",
        rho_b=u.Q_(1.33, "g/cm**3"),
        porosity=u.Q_(0.47),
        theta_s=u.Q_(0.43),
        theta_r=u.Q_(0.078),
        K_sat=u.Q_(2.89e-4, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.036, "1/cm"), n=u.Q_(1.56)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=u.Q_(35.0, "cm")
        ),
    ),
    "Silva et al. (2020) - Loamy sand": soil.Soil(
        name="Silva et al. (2020) - Loamy sand",
        rho_b=u.Q_(1.65, "g/cm**3"),
        porosity=u.Q_(0.44),
        theta_s=u.Q_(0.41),
        theta_r=u.Q_(0.057),
        K_sat=u.Q_(1.23e-3, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0124, "1/cm"), n=u.Q_(2.28)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=u.Q_(35.0, "cm")
        ),
    ),
}

spa_matrix: Dict[soil.Soil, Dict[pfas.PFAS, spa_.SPAParameters]] = {
    soils["Accusand"]: {
        PFASs[
            "PFPeA"
        ]: spa_.KineticFreundlichSorption(  # According to Guo et al. (2022), Table 1
            pfas=PFASs["PFPeA"],
            soil=soils["Accusand"],
            Freundlich_K=u.Q_(0.0211, "(mg/kg)/(mg/L)**0.87"),
            Freundlich_N=u.Q_(0.87),
            frac_instant_adsorption=u.Q_(0.4),
            kinetic_adsorption_rate=u.Q_(5.9, "1/h"),
        ),
        PFASs[
            "PFHxS"
        ]: spa_.KineticFreundlichSorption(  # According to Guo et al. (2022), Table 1
            pfas=PFASs["PFHxS"],
            soil=soils["Accusand"],
            Freundlich_K=u.Q_(0.0213, "(mg/kg)/(mg/L)**0.81"),
            Freundlich_N=u.Q_(0.81),
            frac_instant_adsorption=u.Q_(0.1),
            kinetic_adsorption_rate=u.Q_(3.1, "1/h"),
        ),
        PFASs[
            "PFOA"
        ]: spa_.KineticFreundlichSorption(  # According to Guo et al. (2022), Table 1
            pfas=PFASs["PFOA"],
            soil=soils["Accusand"],
            Freundlich_K=u.Q_(0.1, "mg/kg/(mg/L)**0.87"),
            Freundlich_N=u.Q_(0.87),
            frac_instant_adsorption=u.Q_(0.4),
            kinetic_adsorption_rate=u.Q_(5.9, "1/h"),
        ),
        PFASs[
            "PFOS"
        ]: spa_.KineticFreundlichSorption(  # According to Guo et al. (2022), Table 1
            pfas=PFASs["PFOS"],
            soil=soils["Accusand"],
            Freundlich_K=u.Q_(0.15, "mg/kg/(mg/L)**0.81"),
            Freundlich_N=u.Q_(0.81),
            frac_instant_adsorption=u.Q_(0.1),
            kinetic_adsorption_rate=u.Q_(3.1, "1/h"),
        ),
    },
    soils["Vinton soil"]: {
        PFASs[
            "PFPeA"
        ]: spa_.KineticFreundlichSorption(  # According to Guo et al. (2022), Table 1
            pfas=PFASs["PFPeA"],
            soil=soils["Vinton soil"],
            Freundlich_K=u.Q_(0.122, "(mg/kg)/(mg/L)**0.87"),
            Freundlich_N=u.Q_(0.87),
            frac_instant_adsorption=u.Q_(0.16),
            kinetic_adsorption_rate=u.Q_(0.9, "1/h"),
        ),
        PFASs[
            "PFHxS"
        ]: spa_.KineticFreundlichSorption(  # According to Guo et al. (2022), Table 1
            pfas=PFASs["PFHxS"],
            soil=soils["Vinton soil"],
            Freundlich_K=u.Q_(0.156, "(mg/kg)/(mg/L)**0.77"),
            Freundlich_N=u.Q_(0.77),
            frac_instant_adsorption=u.Q_(0.16),
            kinetic_adsorption_rate=u.Q_(0.9, "1/h"),
        ),
        PFASs[
            "PFOA"
        ]: spa_.KineticFreundlichSorption(  # According to Guo et al. (2022), Table 1
            pfas=PFASs["PFOA"],
            soil=soils["Vinton soil"],
            Freundlich_K=u.Q_(0.58, "(mg/kg)/(mg/L)**0.87"),
            Freundlich_N=u.Q_(0.87),
            frac_instant_adsorption=u.Q_(0.16),
            kinetic_adsorption_rate=u.Q_(0.9, "1/h"),
        ),
        PFASs[
            "PFOS"
        ]: spa_.KineticFreundlichSorption(  # According to Guo et al. (2022), Table 1
            pfas=PFASs["PFOS"],
            soil=soils["Vinton soil"],
            Freundlich_K=u.Q_(1.11, "(mg/kg)/(mg/L)**0.77"),
            Freundlich_N=u.Q_(0.77),
            frac_instant_adsorption=u.Q_(0.16),
            kinetic_adsorption_rate=u.Q_(0.9, "1/h"),
        ),
    },
    soils["Schoonenburgse Heuvel - sand"]: {
        PFASs[
            "PFOA"
        ]: spa_.LinearSorption(  # According to de Jong, 2022, source unknown
            pfas=PFASs["PFOA"],
            soil=soils["Schoonenburgse Heuvel - sand"],
            Kd_=u.Q_(1.083, "cm^3/g"),
        ),
        PFASs[
            "PFOS"
        ]: spa_.LinearSorption(  # According to de Jong, 2022, source unknown
            pfas=PFASs["PFOS"],
            soil=soils["Schoonenburgse Heuvel - sand"],
            Kd_=u.Q_(2.067, "cm^3/g"),
        ),
    },
}

PFASRegistry = r.PFASRegistry(
    pfas=PFASs,
    soil=soils,
    spa=spa_matrix,
)
