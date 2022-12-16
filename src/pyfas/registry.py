import dataclasses
from typing import Dict, Generator, Tuple
import warnings

from . import pfas as p
from . import soil as s
from . import solid_phase_adsorption as spa_
from . import air_water_adsorption as awa_


@dataclasses.dataclass
class PFASRegistry:
    pfas: Dict[str, p.PFAS]
    soil: Dict[str, s.Soil]
    spa: Dict[s.Soil, Dict[p.PFAS, spa_.SPAParameters]]
    awa: Dict[p.PFAS, awa_.AWAParameters]

    def get_parameters(
        self, pfas_name: str, soil_name: str
    ) -> Tuple[p.PFAS, s.Soil, spa_.SPAParameters, awa_.AWAParameters]:
        try:
            pfas = self.pfas[pfas_name]
        except KeyError as e:
            raise ValueError(f"PFAS {pfas_name} not found in registry.") from e
        try:
            soil = self.soil[soil_name]
        except KeyError as e:
            raise ValueError(f"Soil {soil_name} not found in registry.") from e
        try:
            spa = self.spa[soil][pfas]
        except KeyError as e:
            warnings.warn(
                f"PFAS {pfas_name} is not registered for soil {soil_name}, returning a distributed sorption model."
            )
            spa = spa_.FabregatPalauSorption(
                pfas,
                soil,
            )
        try:
            awa = self.awa[pfas]
        except KeyError as e:
            warnings.warn(
                f"PFAS {pfas_name} is not registered for air-water partitioning, returning a distributed air-water partitioning model."
            )
            awa = awa_.LeSorption.from_pfas(pfas)
        return pfas, soil, spa, awa

    def __iter__(
        self,
    ) -> Generator[
        Tuple[p.PFAS, s.Soil, spa_.SPAParameters, awa_.AWAParameters], None, None
    ]:
        for pfas in self.pfas.values():
            for soil in self.soil.values():
                if soil in self.spa:
                    if pfas in self.spa[soil] and pfas in self.awa:
                        yield pfas, soil, self.spa[soil][pfas], self.awa[pfas]
