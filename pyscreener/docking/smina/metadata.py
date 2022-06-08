from dataclasses import dataclass
from os import PathLike
import shlex
from typing import Iterable, Optional, Union

from pyscreener.docking.metadata import SimulationMetadata

@dataclass(repr=True, eq=False)
class SminaMetadata(SimulationMetadata):
    """

    Attributes
    ---------
    exhaustiveness: int
        the exhaustiveness of the global search. Larger values are more exhaustive
    num_modes: int
        the number of output modes
    energy_range: float
        the maximum energy difference (in kcal/mol) between the best and worst output binding modes
    extra : List[str]
        additional arguments that will be passed to the docking calculation
    prepared_ligand: Optional[Path]
    prepared_receptor: Optional[Path]

    Parmeters
    ---------
    exhaustiveness: int, default=8
    num_modes: int, default=9
    energy_range: float, default=3.
    extra : str, default=""
        a string containing the additional command line arguments to pass to a run of a vina-type
        software for options not contained within the default metadata. E.g. for a run of Smina, extra="--force_cap ARG" or for PSOVina, extra="-w ARG"
    prepared_ligand: Optional[Union[str, Path]] = None,
    prepared_receptor: Optional[Union[str, Path]] = None
    """

    exhaustiveness: int = 8
    num_modes: int = 9
    energy_range: float = 3.0
    extra: Union[str, Iterable[str]] = ""
    prepared_ligand: Optional[PathLike] = None
    prepared_receptor: Optional[PathLike] = None

    def __post_init__(self):
        self.extra = shlex.split(self.extra) if isinstance(self.extra, str) else self.extra
