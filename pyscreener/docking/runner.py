from abc import ABC, abstractmethod
from typing import Optional, Sequence

from pyscreener.docking.result import Result
from pyscreener.docking.sim import Simulation
from pyscreener.docking.metadata import SimulationMetadata


class DockingRunner(ABC):
    @classmethod
    @property
    @abstractmethod
    def is_multithreaded(cls) -> bool:
        """Is this docking program able to leverage multiple CPU cores?"""

    @staticmethod
    @abstractmethod
    def prepare_receptor(sim: Simulation) -> Simulation:
        """Prepare the receptor file(s) for the given simulation"""

    @staticmethod
    @abstractmethod
    def prepare_ligand(sim: Simulation) -> bool:
        """Prepare the ligand file(s) for the given simulation and return True upon success"""

    @staticmethod
    @abstractmethod
    def run(sim: Simulation) -> Optional[Sequence[float]]:
        """Run the given simulation and return the score(s) of the docked conformers"""

    @staticmethod
    @abstractmethod
    def prepare_and_run(sim: Simulation) -> Optional[Result]:
        """Prepare the ligand file then run the given simulation. Roughly equivlaent to `prepare_ligand()` followed by `run()` but returns the Result object for the Simulation
        rather than the scores of the conformers"""

    @staticmethod
    def validate_metadata(md: SimulationMetadata):
        """Validate the metadata of the simulation. E.g., ensure that the specified software is
        installed for Vina-type screens."""
