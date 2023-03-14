from typing import Dict

from . import _typing as t
from . import air_water_adsorption as awa_
from . import pfas
from . import registry as r
from . import soil
from . import solid_phase_adsorption as spa_
from . import units as u

PFASs = {
    "TFA": pfas.PFAS(
        name="TFA",
        M=u.Q_(114.02, "g/mol"),
        diffusivity=pfas.estimate_diffusion_coefficient(1),
        n_CFx=u.Q_(1),
        n_COO=u.Q_(1),
    ),
    "PFBA": pfas.PFAS(
        name="PFBA",
        M=u.Q_(214.0, "g/mol"),
        diffusivity=pfas.estimate_diffusion_coefficient(3),
        K_oc=u.Q_(2.9, "L/kg"),
        K_sc=u.Q_(0.43, "L/kg"),
        n_CFx=u.Q_(3),
        n_COO=u.Q_(1),
    ),
    "PFBS": pfas.PFAS(
        name="PFBS",
        M=u.Q_(300.1, "g/mol"),  # According to Guo et al. (2022), Table 2
        diffusivity=t.Diffusivity(
            11e-6, "cm**2/s"
        ),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(11.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(0.44, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CFx=u.Q_(4),
        n_SO3=u.Q_(1),
    ),
    "PFPeA": pfas.PFAS(
        name="PFPeA",
        M=u.Q_(264.05, "g/mol"),  # According to Guo et al. (2022), Table 2
        diffusivity=t.Diffusivity(
            12e-6, "cm**2/s"
        ),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(15.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(0.46, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CFx=u.Q_(4),
        n_COO=u.Q_(1),
    ),
    "PFHxA": pfas.PFAS(
        name="PFHxA",
        M=u.Q_(314.05, "g/mol"),  # According to Guo et al. (2022), Table 2
        diffusivity=t.Diffusivity(
            7.8e-6, "cm**2/s"
        ),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(15.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(0.46, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CFx=u.Q_(5),
        n_COO=u.Q_(1),
    ),
    "PFHxS": pfas.PFAS(
        name="PFHxS",
        M=u.Q_(400.12, "g/mol"),  # According to Guo et al. (2022), Table 2
        diffusivity=t.Diffusivity(
            4.5e-6, "cm**2/s"
        ),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(50.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(1.2, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CFx=u.Q_(6),
        n_SO3=u.Q_(1),
    ),
    "PFHpA": pfas.PFAS(
        name="PFHpA",
        M=u.Q_(364.06, "g/mol"),  # According to Guo et al. (2022), Table 2
        diffusivity=t.Diffusivity(
            9.3e-6, "cm**2/s"
        ),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(50.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=None,
        K_mo=None,
        n_CFx=u.Q_(6),
        n_COO=u.Q_(1),
    ),
    "PFOA": pfas.PFAS(
        name="PFOA",
        M=u.Q_(414.07, "g/mol"),  # According to Guo et al. (2022), Table 2
        diffusivity=t.Diffusivity(
            4.9e-6, "cm**2/s"
        ),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(107.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(3.3, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CFx=u.Q_(7),
        n_COO=u.Q_(1),
    ),
    "PFOS": pfas.PFAS(
        name="PFOS",
        M=u.Q_(500.13, "g/mol"),  # According to Guo et al. (2022), Table 2
        diffusivity=t.Diffusivity(
            5.4e-6, "cm**2/s"
        ),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(609.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(9.4, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CFx=u.Q_(8),
        n_SO3=u.Q_(1),
    ),
    "PFNA": pfas.PFAS(
        name="PFNA",
        M=u.Q_(464.08, "g/mol"),  # According to Guo et al. (2022), Table 2
        diffusivity=t.Diffusivity(
            2.93e-6, "cm**2/s"
        ),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(324.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(2.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CFx=u.Q_(8),
        n_COO=u.Q_(1),
    ),
    "PFDA": pfas.PFAS(
        name="PFDA",
        M=u.Q_(514.08, "g/mol"),  # According to Guo et al. (2022), Table 2
        diffusivity=t.Diffusivity(
            2.27e-6, "cm**2/s"
        ),  # According to Guo et al. (2022), Table 2
        K_oc=u.Q_(604.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_sc=u.Q_(14.0, "L/kg"),  # According to Fabregat-Palau et al. (2021), Table 3
        K_mo=None,
        n_CFx=u.Q_(9),
        n_COO=u.Q_(1),
    ),
    "HFPO-DA": pfas.PFAS(
        name="HFPO-DA",
        M=u.Q_(
            330.05, "g/mol"
        ),  # According to https://pubchem.ncbi.nlm.nih.gov/compound/114481
        diffusivity=pfas.estimate_diffusion_coefficient(5),
        n_CFx=u.Q_(5),
        n_COO=u.Q_(1),
        n__O_=u.Q_(1),
    ),
}
soils = {
    "Accusand": soil.Soil(  # According to Guo et al. (2020), Table 1
        name="Accusand",
        rho_b=soil.BulkDensity(1.65, "g/cm**3"),
        porosity=soil.Porosity(0.294),
        theta_s=soil.SaturatedWaterContent(0.294),
        theta_r=soil.ResidualWaterContent(0.015),
        K_sat=soil.SaturatedHydraulicConductivity(2.0964e-2, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.04479, "1/cm"), n=u.Q_(4.0)
        ),
        tracer_fit=soil.TracerFitParameters(
            x0=u.Q_(633.96, "cm**2/cm**3"),
            x1=u.Q_(-1182.5, "cm**2/cm**3"),
            x2=u.Q_(548.54, "cm**2/cm**3"),
        ),
        soil_roughness_multiplier=u.Q_(4.15),
        f_oc=soil.FractionOrganicCarbon(
            0.04, "percent"
        ),  # According to Guo et al. (2020), Section 4
        f_mo=soil.FractionMetalOxides(
            (14.0 + 2.5 + 12.0), "ug/g"  # FeOx + MnOx + AlOx
        ),  # According to Guo et al. (2020), Section 4
        f_clay=soil.FractionClay(
            0.0, "percent"
        ),  # According to Guo et al. (2020), Section 4
        f_silt=soil.FractionSilt(
            0.0, "percent"
        ),  # According to Guo et al. (2020), Section 4
    ),
    "Vinton soil": soil.Soil(  # According to Guo et al. (2020), Table 1
        name="Vinton soil",
        rho_b=soil.BulkDensity(1.627, "g/cm**3"),
        porosity=soil.Porosity(0.395),
        theta_s=soil.SaturatedWaterContent(0.395),
        theta_r=soil.ResidualWaterContent(0.056),
        K_sat=soil.SaturatedHydraulicConductivity(1.17e-3, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.02178, "1/cm"), n=u.Q_(3.451)
        ),
        tracer_fit=soil.TracerFitParameters(
            x0=u.Q_(1543.6, "cm**2/cm**3"),
            x1=u.Q_(-2848.6, "cm**2/cm**3"),
            x2=u.Q_(1305.0, "cm**2/cm**3"),
        ),
        soil_roughness_multiplier=u.Q_(4.15),
        f_oc=soil.FractionOrganicCarbon(
            0.1, "percent"
        ),  # According to Guo et al. (2020), Section 4
        f_mo=soil.FractionMetalOxides(
            0.0, "percent"
        ),  # According to Guo et al. (2020), Section 4
        f_clay=soil.FractionClay(
            4.7, "percent"
        ),  # According to Guo et al. (2020), Section 4
        f_silt=soil.FractionSilt(
            0.0, "percent"
        ),  # According to Guo et al. (2020), Section 4
    ),
    "Schoonenburgse Heuvel - sand": soil.Soil(  # According to De Jong (unpublished), Appendix A
        name="Schoonenburgse Heuvel - sand",
        rho_b=soil.BulkDensity(1.5, "g/cm**3"),
        porosity=soil.Porosity(0.427),
        theta_s=soil.SaturatedWaterContent(0.427),
        theta_r=soil.ResidualWaterContent(0.02),
        K_sat=soil.SaturatedHydraulicConductivity(3.61e-4, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0217, "1/cm"), n=u.Q_(1.735)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=soil.Dispersivity(5.0, "cm")
        ),
    ),
    "Schoonenburgse Heuvel - peat": soil.Soil(  # According to De Jong (unpublished), Appendix A
        name="Schoonenburgse Heuvel - peat",
        rho_b=soil.BulkDensity(0.23, "g/cm**3"),
        porosity=soil.Porosity(0.85),
        theta_s=soil.SaturatedWaterContent(0.849),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(3.93519e-05, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0119, "1/cm"), n=u.Q_(1.272)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=soil.Dispersivity(10.0, "cm")
        ),
    ),
    "Silva et al. (2020) - loam": soil.Soil(  # According to De Jong (unpublished), Appendix A
        name="Silva et al. (2020) - loam",
        rho_b=soil.BulkDensity(1.33, "g/cm**3"),
        porosity=soil.Porosity(0.47),
        theta_s=soil.SaturatedWaterContent(0.43),
        theta_r=soil.ResidualWaterContent(0.078),
        K_sat=soil.SaturatedHydraulicConductivity(2.89e-4, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.036, "1/cm"), n=u.Q_(1.56)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=soil.Dispersivity(35.0, "cm")
        ),
    ),
    "Silva et al. (2020) - Loamy sand": soil.Soil(  # According to De Jong (unpublished), Appendix A
        name="Silva et al. (2020) - Loamy sand",
        rho_b=soil.BulkDensity(1.65, "g/cm**3"),
        porosity=soil.Porosity(0.44),
        theta_s=soil.SaturatedWaterContent(0.41),
        theta_r=soil.ResidualWaterContent(0.057),
        K_sat=soil.SaturatedHydraulicConductivity(1.23e-3, "cm/s"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0124, "1/cm"), n=u.Q_(2.28)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        dispersivity=soil.longitudinal_dispersivity_constant(
            dispersivity=soil.Dispersivity(35.0, "cm")
        ),
    ),
    "Staring-B01": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-B01",
        rho_b=soil.BulkDensity(1.5, "g/cm**3"),  # Unknown
        porosity=soil.Porosity(0.427),
        theta_s=soil.SaturatedWaterContent(0.427),
        theta_r=soil.ResidualWaterContent(0.02),
        K_sat=soil.SaturatedHydraulicConductivity(31.23, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0217, "1/cm"), n=u.Q_(1.735), l=u.Q_(0.981)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(5.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(7.5, "percent"),
    ),
    "Staring-O01": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O01",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.366),
        ),
        porosity=soil.Porosity(0.366),
        theta_s=soil.SaturatedWaterContent(0.366),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(22.32, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0160, "1/cm"), n=u.Q_(2.163), l=u.Q_(2.868)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(5.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O02": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O02",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.387),
        ),
        porosity=soil.Porosity(0.387),
        theta_s=soil.SaturatedWaterContent(0.387),
        theta_r=soil.ResidualWaterContent(0.02),
        K_sat=soil.SaturatedHydraulicConductivity(22.76, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0161, "1/cm"), n=u.Q_(1.524), l=u.Q_(2.440)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(14.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O03": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O03",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.340),
        ),
        porosity=soil.Porosity(0.340),
        theta_s=soil.SaturatedWaterContent(0.340),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(12.37, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0172, "1/cm"), n=u.Q_(1.703), l=u.Q_(0.0)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(25.5, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O04": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O04",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.364),
        ),
        porosity=soil.Porosity(0.364),
        theta_s=soil.SaturatedWaterContent(0.364),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(25.81, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0136, "1/cm"), n=u.Q_(1.488), l=u.Q_(2.179)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(41.5, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O05": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O05",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.337),
        ),
        porosity=soil.Porosity(0.337),
        theta_s=soil.SaturatedWaterContent(0.337),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(17.42, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0303, "1/cm"), n=u.Q_(2.888), l=u.Q_(0.074)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O06": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O06",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.333),
        ),
        porosity=soil.Porosity(0.333),
        theta_s=soil.SaturatedWaterContent(0.333),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(32.83, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0160, "1/cm"), n=u.Q_(1.289), l=u.Q_(-1.010)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(25.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O07": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O07",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.513),
        ),
        porosity=soil.Porosity(0.513),
        theta_s=soil.SaturatedWaterContent(0.513),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(37.55, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0120, "1/cm"), n=u.Q_(1.153), l=u.Q_(-2.013)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(41.5, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O08": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O08",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(10.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.454),
        ),
        porosity=soil.Porosity(0.454),
        theta_s=soil.SaturatedWaterContent(0.454),
        theta_r=soil.ResidualWaterContent(0.0),
        K_sat=soil.SaturatedHydraulicConductivity(8.64, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0113, "1/cm"), n=u.Q_(1.346), l=u.Q_(-0.904)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(10.0, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O09": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O09",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(15.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.458),
        ),
        porosity=soil.Porosity(0.458),
        theta_s=soil.SaturatedWaterContent(0.458),
        theta_r=soil.ResidualWaterContent(0.0),
        K_sat=soil.SaturatedHydraulicConductivity(3.77, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0097, "1/cm"), n=u.Q_(1.376), l=u.Q_(-1.013)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(15.0, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O10": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O10",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(21.5, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.472),
        ),
        porosity=soil.Porosity(0.472),
        theta_s=soil.SaturatedWaterContent(0.472),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(2.30, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0100, "1/cm"), n=u.Q_(1.246), l=u.Q_(-0.793)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(21.5, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O11": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O11",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(30.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.444),
        ),
        porosity=soil.Porosity(0.444),
        theta_s=soil.SaturatedWaterContent(0.444),
        theta_r=soil.ResidualWaterContent(0.0),
        K_sat=soil.SaturatedHydraulicConductivity(2.12, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0143, "1/cm"), n=u.Q_(1.126), l=u.Q_(2.357)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(30.0, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O12": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O12",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(42.5, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.561),
        ),
        porosity=soil.Porosity(0.561),
        theta_s=soil.SaturatedWaterContent(0.561),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(1.08, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0088, "1/cm"), n=u.Q_(1.158), l=u.Q_(-3.172)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(42.5, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O13": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O13",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(75.0, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.573),
        ),
        porosity=soil.Porosity(0.573),
        theta_s=soil.SaturatedWaterContent(0.573),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(9.69, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0279, "1/cm"), n=u.Q_(1.080), l=u.Q_(-6.091)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(75.0, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O14": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O14",
        rho_b=soil.bulk_density_Poelman1974(
            f_clay=soil.FractionClay(0.0, "percent"),
            f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
            porosity=soil.Porosity(0.394),
        ),
        porosity=soil.Porosity(0.394),
        theta_s=soil.SaturatedWaterContent(0.394),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(2.50, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0033, "1/cm"), n=u.Q_(1.617), l=u.Q_(-0.514)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(67.5, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O15": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O15",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(92.5, "percent"),
            soil.FractionOrganicCarbon(1.5, "percent"),
            soil.Porosity(0.410),
        ),
        porosity=soil.Porosity(0.410),
        theta_s=soil.SaturatedWaterContent(0.410),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(2.79, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0078, "1/cm"), n=u.Q_(1.287), l=u.Q_(0.000)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(92.5, "percent"),
        f_oc=soil.FractionOrganicCarbon(1.5, "percent"),
    ),
    "Staring-O16": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O16",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(67.5, "percent"),
            soil.Porosity(0.889),
        ),
        porosity=soil.Porosity(0.889),
        theta_s=soil.SaturatedWaterContent(0.889),
        theta_r=soil.ResidualWaterContent(0.0),
        K_sat=soil.SaturatedHydraulicConductivity(1.46, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0097, "1/cm"), n=u.Q_(1.364), l=u.Q_(-0.665)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(67.5, "percent"),
    ),
    "Staring-O17": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O17",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(67.5, "percent"),
            soil.Porosity(0.849),
        ),
        porosity=soil.Porosity(0.849),
        theta_s=soil.SaturatedWaterContent(0.849),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(3.40, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0119, "1/cm"), n=u.Q_(1.272), l=u.Q_(-1.249)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(67.5, "percent"),
    ),
    "Staring-O18": soil.Soil(  # According to Heinen et al. (2020), Tables 2 and 3
        name="Staring-O18",
        rho_b=soil.bulk_density_Poelman1974(
            soil.FractionClay(0.0, "percent"),
            soil.FractionOrganicCarbon(25.0, "percent"),
            soil.Porosity(0.580),
        ),
        porosity=soil.Porosity(0.580),
        theta_s=soil.SaturatedWaterContent(0.580),
        theta_r=soil.ResidualWaterContent(0.01),
        K_sat=soil.SaturatedHydraulicConductivity(35.97, "cm/day"),
        van_genuchten=soil.VanGenuchtenParameters(
            alpha=u.Q_(0.0127, "1/cm"), n=u.Q_(1.316), l=u.Q_(-0.786)
        ),
        tracer_fit=None,
        soil_roughness_multiplier=u.Q_(4.15),
        f_clay=soil.FractionClay(0.0, "percent"),
        f_silt=soil.FractionSilt(0.0, "percent"),
        f_oc=soil.FractionOrganicCarbon(25.0, "percent"),
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
awa_matrix: Dict[pfas.PFAS, awa_.AWAParameters] = {
    PFASs["PFPeA"]: awa_.SzyszkowskiSorption(
        a=(u.Q_(3168.6, "mg/L") / u.Q_(264.05, "g/mol")).to(
            "µmol/L"
        ),  # According to Guo et al. (2022), Table 2
        b=u.Q_(0.22),  # According to Guo et al. (2022), Table 2
    ),
    PFASs["PFHxA"]: awa_.SzyszkowskiSorption(
        a=(u.Q_(1350.42, "mg/L") / u.Q_(314.05, "g/mol")).to(
            "µmol/L"
        ),  # According to Guo et al. (2022), Table 2
        b=u.Q_(0.21),  # According to Guo et al. (2022), Table 2
    ),
    PFASs["PFHpA"]: awa_.SzyszkowskiSorption(
        a=(u.Q_(345.86, "mg/L") / u.Q_(364.06, "g/mol")).to(
            "µmol/L"
        ),  # According to Guo et al. (2022), Table 2
        b=u.Q_(0.22),  # According to Guo et al. (2022), Table 2
    ),
    PFASs["PFOA"]: awa_.SzyszkowskiSorption(
        a=u.Q_(62.11, "mg/L")
        / u.Q_(414.07, "g/mol"),  # According to Guo et al. (2022), Table 2
        b=u.Q_(0.19),  # According to Guo et al. (2022), Table 2
    ),
    PFASs["PFNA"]: awa_.SzyszkowskiSorption(
        a=u.Q_(5.11, "mg/L")
        / u.Q_(464.08, "g/mol"),  # According to Guo et al. (2022), Table 2
        b=u.Q_(0.16),  # According to Guo et al. (2022), Table 2
    ),
    PFASs["PFDA"]: awa_.SzyszkowskiSorption(
        a=u.Q_(3.7, "mg/L")
        / u.Q_(514.08, "g/mol"),  # According to Guo et al. (2022), Table 2
        b=u.Q_(0.17),  # According to Guo et al. (2022), Table 2
    ),
    PFASs["PFBS"]: awa_.SzyszkowskiSorption(
        a=u.Q_(2400.8, "mg/L")
        / u.Q_(300.1, "g/mol"),  # According to Guo et al. (2022), Table 2
        b=u.Q_(0.15),  # According to Guo et al. (2022), Table 2
    ),
    PFASs["PFHxS"]: awa_.SzyszkowskiSorption(
        a=u.Q_(160.05, "mg/L")
        / u.Q_(400.12, "g/mol"),  # According to Guo et al. (2022), Table 2
        b=u.Q_(0.14),  # According to Guo et al. (2022), Table 2
    ),
    PFASs["PFOS"]: awa_.SzyszkowskiSorption(
        a=u.Q_(3.65, "mg/L")
        / u.Q_(500.13, "g/mol"),  # According to Guo et al. (2022), Table 2
        b=u.Q_(0.12),  # According to Guo et al. (2022), Table 2
    ),
}

PFASRegistry = r.PFASRegistry(
    pfas=PFASs,
    soil=soils,
    spa=spa_matrix,
    awa=awa_matrix,
)
