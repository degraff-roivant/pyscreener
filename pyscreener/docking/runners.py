from __future__ import annotations

from abc import ABC, abstractmethod
from os import PathLike
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
    def check_environment(md: SimulationMetadata):
        """Check that the environment is properly set up for the given runner and metadata.
        
        Raises
        ------
        InvalidEnvironmentError
            if the environment is not set up properly
        """


class BatchDockingRunner(DockingRunner):
    @staticmethod
    @abstractmethod
    def batch_prepare_and_run(sims: Sequence[Simulation]) -> list[Optional[Result]]:
        pass

    @staticmethod
    @abstractmethod
    def batch_prepare_ligand(sims: Sequence[Simulation]) -> bool:
        pass

    @staticmethod
    @abstractmethod
    def batch_run(sims: Sequence[Simulation]) -> Optional[list[float]]:
        pass

    @staticmethod
    @abstractmethod
    def batch_parse_logfile(logfile: PathLike):
        pass
