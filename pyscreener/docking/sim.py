from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

from pyscreener.exceptions import InvalidResultError, NotSimulatedError
from pyscreener.utils import Reduction
from pyscreener.docking.metadata import SimulationMetadata
from pyscreener.docking.result import Result


@dataclass(repr=True, eq=False)
class Simulation:
    """

    Attributes
    ---------
    smi : str
        the SMILES string of the ligand that will be docked
    receptor : Path
        the path of the receptor to use for docking
    metadata : CalculationMetadata
        the parameters with which to prepare and run the simulation
    name : str
        the name to use when creating the ligand input file and output files
    input_file : Optional[Path]
        the filepath of an arbitrary molecular input file containing a single molecule
    in_path : Union[str, Path]
        the path under which input will be placed
    out_path: Union[str, Path]
        the path under which output will be placed
    reduction : str
        the mode used to calculate a score for an individual docking calculation given multiple
        output scored conformations
    k : int
        the number of top scores to use if calculating an average
    result : Optional[Mapping]
        the result of the docking calculation. None if the calculation has not been performed yet.

    Parmeters
    ---------
    smi : str
    receptor : PathLike
    metadata : CalculationMetadata
    name : str, default='ligand'
    input_file : Optional[str], default=None
    in_path : PathLike, default='.'
    out_path : PathLike, default='.'
    reduction : str, default=ScoreMode.BEST
    k : int, default=1
    result : Optional[Result], default=None
    """

    smi: str
    receptor: Path
    center: Tuple[float, float, float]
    size: Tuple[float, float, float]
    metadata: SimulationMetadata
    ncpu: int = 1
    name: str = "ligand"
    input_file: Optional[Path] = None
    in_path: Path = Path(".")
    out_path: Path = Path(".")
    reduction: Reduction = Reduction.BEST
    k: int = 1
    result: Optional[Result] = None

    def __post_init__(self):
        self.receptor = Path(self.receptor).absolute()
        if self.input_file is not None:
            self.input_file = Path(self.input_file)
        self.in_path = Path(self.in_path)
        self.out_path = Path(self.out_path)

    @property
    def score(self) -> Optional[float]:
        """the docking score of this calculation

        Raises
        ------
        NotSimulatedError
            if this calculation has not been run yet
        """
        if self.result is None:
            raise NotSimulatedError("Simulation has not been run!")

        try:
            return self.result.score
        except AttributeError:
            raise InvalidResultError('No attribute: "score" in Result object (self.result)!')
