from __future__ import annotations

from itertools import takewhile
from os import PathLike
from pathlib import Path
import re
import shutil
import subprocess as sp
from typing import Optional, Sequence
import warnings

import numpy as np
from rdkit.Chem import AllChem as Chem
import ray

from pyscreener import utils
from pyscreener.docking.utils import parse_optional_float
from pyscreener.exceptions import MissingExecutableError
from pyscreener.warnings import SimulationFailureWarning
from pyscreener.docking.sim import Simulation
from pyscreener.docking.runners import BatchDockingRunner
from pyscreener.docking.result import Result
from pyscreener.docking.smina.metadata import SminaMetadata


class SminaRunner(BatchDockingRunner):
    @classmethod
    @property
    def is_multithreaded(cls) -> bool:
        return True

    @staticmethod
    def prepare(sim: Simulation) -> Simulation:
        SminaRunner.prepare_receptor(sim)
        SminaRunner.prepare_ligand(sim)

        return sim

    @staticmethod
    def prepare_receptor(sim: Simulation) -> Simulation:
        """Set the `prepared_receptor` attribute of the metadata appropriately"""
        dest = Path(sim.in_path) / sim.receptor.name

        shutil.copy(str(sim.receptor), str(dest))
        sim.metadata.prepared_receptor = dest

        return sim

    @staticmethod
    def prepare_and_run(sim: Simulation) -> Optional[Result]:
        if not SminaRunner.prepare_ligand(sim):
            return None

        _ = SminaRunner.run(sim)

        return sim.result

    @staticmethod
    def prepare_ligand(sim: Simulation) -> bool:
        if sim.smi is not None:
            return SminaRunner.prepare_from_smi(sim)

        return SminaRunner.prepare_from_file(sim)

    @staticmethod
    def prepare_from_smi(sim: Simulation) -> bool:
        """Prepare the ligand PDQBT file from its SMILES string

        If successful, sets the `prepared_ligand` attribute of the metadata  and return True.
        Otherwise, do nothing and return False.

        Parameters
        ----------
        sim : Simulation

        Returns
        -------
        bool
            whether the ligand preparation succeeded
        """
        sdf = Path(sim.in_path) / f"{sim.name}.sdf"
        writer = Chem.SDWriter(str(sdf))

        mol = Chem.MolFromSmiles(sim.smi)
        if mol is None:
            return False

        writer.write(mol)
        sim.metadata.prepared_ligand = sdf

        return True

    @staticmethod
    def prepare_from_file(sim: Simulation) -> bool:
        """Prepare the ligand input file from its input chemical file, retaining the input
        geometry and sets the `prepared_ligand` attribute of the metadata in the process. If ligand
        preparation fails for any reason, do nothing and return False.

        NOTE: the input chemical file must be a file that is readable by Smina

        Parameters
        ----------
        sim : Simulation

        Returns
        -------
        bool
            whether the ligand preparation succeeded
        """
        sim.metadata.prepared_ligand = sim.input_file

        return True

    @staticmethod
    def run(sim: Simulation) -> Optional[np.ndarray]:
        """run the given ligand using the specified vina-type docking program and parameters

        If the simulation is not possible due to the prepared inputs not being set beforehand, do
        nothing and return None. Otherwise, if the simulation itself fails, set the `result`
        attribute of the simulation with a score of `None` **and** return `None`.

        Returns
        -------
        scores : Optional[list[float]]
            the conformer scores parsed from the log file. None if no scores were parseable from
            the logfile due to simulation failure.
        """
        if sim.metadata.prepared_receptor is None or sim.metadata.prepared_ligand is None:
            return None

        p_ligand = Path(sim.metadata.prepared_ligand)
        ligand_name = p_ligand.stem

        name = f"{Path(sim.receptor).stem}_{ligand_name}"

        argv, _, log = SminaRunner.build_argv(
            sim.metadata.prepared_ligand,
            sim.metadata.prepared_receptor,
            sim.center,
            sim.size,
            sim.ncpu,
            sim.metadata.exhaustiveness,
            sim.metadata.num_modes,
            sim.metadata.energy_range,
            name,
            Path(sim.out_path),
            sim.metadata.extra,
        )

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            warnings.warn(f'Message: {ret.stderr.decode("utf-8")}', SimulationFailureWarning)

        scores = SminaRunner.parse_logfile(log)
        if scores is None:
            score = None
        else:
            score = utils.reduce_scores(np.array(scores), sim.reduction, k=sim.k)

        sim.result = Result(sim.smi, name, re.sub("[:,.]", "", ray.state.current_node_id()), score)

        return scores

    @staticmethod
    def build_argv(
        ligand: str,
        receptor: str,
        center: tuple[float, float, float],
        size: tuple[float, float, float] = (10, 10, 10),
        ncpu: int = 1,
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0,
        autobox_ligand: Optional[PathLike] = None,
        autobox_add: float = 4.,
        name: Optional[str] = None,
        path: Path = Path("."),
        extra: Optional[list[str]] = None,
    ) -> tuple[list[str], Path, Path]:
        """Builds the argument vector to run a vina-type docking program

        Parameters
        ----------
        ligand : str
            the filename of the input ligand PDBQT file
        receptor : str
            the filename of the input receptor PDBQT file
        center : tuple[float, float, float]
            the x-, y-, and z-coordinates of the center of the search box
        size : tuple[float, float, float], default=(10, 10, 10)
            the  x-, y-, and z-radii, respectively, of the search box
        ncpu : int, default=1
            the number of cores to allocate to the docking program
        exhaustiveness: int
            the exhaustiveness of the global search. Larger values are more exhaustive
        num_modes: int
            the number of output modes
        energy_range: float
            the maximum energy difference (in kcal/mol) between the best and worst output binding
            modes
        autobox_ligand: Optional[PathLike], default=None
            the ligand to autobox.
            NOTE: ligand autoboxing will take priority over input docking box parameters!
        autobox_add: float, default=4.
            the amount of buffer to add when autoboxing a ligand
        extra : Optional[list[str]]
            additional command line arguments that will be passed to the docking calculation
        name : string, default=<receptor>_<ligand>)
            the base name to use for both the log and out files
        path : Path, default=Path('.')
            the path under which both the log and out files should be written
        extra : Optional[list[str]], default=None
            additional command line arguments to pass to each run

        Returns
        -------
        argv : list[str]
            the argument vector with which to run an instance of a vina-type
            docking program
        out : Path
            the filepath of the out file which the docking program will write to
        log : Path
            the filepath of the log file which the docking program will write to
        """
        name = name or (Path(receptor).stem + "_" + Path(ligand).stem)
        extra = extra or []

        out = path / f"smina_{name}_out.pdb"
        log = path / f"smina_{name}.log"

        argv = [
            "smina",
            f"--receptor={receptor}",
            f"--ligand={ligand}",
            f"--cpu={ncpu}",
            f"--out={out}",
            f"--log={log}",
            f"--exhaustiveness={exhaustiveness}",
            f"--num_modes={num_modes}",
            f"--energy_range={energy_range}",
        ]
        if autobox_ligand is not None:
            argv.extend([f"--autobox_ligand={autobox_ligand}", f"--autobox_add={autobox_add}"])
        else:
            argv.extend([
                f"--center_x={center[0]}",
                f"--center_y={center[1]}",
                f"--center_z={center[2]}",
                f"--size_x={size[0]}",
                f"--size_y={size[1]}",
                f"--size_z={size[2]}"
            ])
        argv.extend(extra)

        return argv, out, log

    @staticmethod
    def parse_logfile(logfile: PathLike) -> Optional[list[float]]:
        """parse a Smina-type log file for the scores of the binding modes

        Parameters
        ----------
        logfile : PathLike
            the path to a Vina-type log file

        Returns
        -------
        Optional[list[float]]
            the scores of the docked binding modes in the ordering of the log file. None if no
            scores were parsed or the log file was unparsable
        """
        SCORE_LINE_PATTERN = r"\d\s+\-?\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+"
        pattern = re.compile(SCORE_LINE_PATTERN)

        try:
            with open(logfile) as fid:
                for line in fid:
                    if pattern.match(line) is not None:
                        break
                score_lines = [line]
                score_lines.extend(
                    list(takewhile(lambda line: pattern.match(line) is not None, fid))
                )
        except OSError:
            return None

        scores = np.array([parse_optional_float(line.split()[1]) for line in score_lines])
        return scores if len(scores) > 0 else None

    @staticmethod
    def parse_outfile(outfile: PathLike) -> Optional[np.ndarray]:
        """parse a Smina-type output file for the scores of the binding modes

        Paramaters
        ----------
        outfile : PathLike
            the filepath a vina-type output file

        Returns
        -------
        Optional[list[float]]
            the scores of the binding in the ordering of the output file. None if no scores were
            parsed or the log file was unparseable
        """
        try:
            with open(outfile) as fid:
                score_lines = [
                    line for line in fid.readlines() if "REMARK minimizedAffinity" in line
                ]
        except OSError:
            return None

        scores = np.array([parse_optional_float(line.split()[3]) for line in score_lines])
        return scores if len(scores) > 0 else None

    @staticmethod
    def check_environment(metadata: SminaMetadata):
        if shutil.which("smina") is None:
            raise MissingExecutableError(
                "Could not find `smina` on PATH! "
                "See https://github.com/coleygroup/pyscreener#adding-an-executable-to-your-path for more information."
            )

    @staticmethod
    def batch_prepare_and_run(sims: Sequence[Simulation]) -> list[Optional[Result]]:
        SminaRunner.batch_prepare_ligand(sims)
        SminaRunner.batch_run(sims)

        return [sim.result for sim in sims]

    @staticmethod
    def batch_prepare_ligand(sims: Sequence[Simulation]) -> bool:
        """Prepare a batch of ligands. NOTE: only works for SMILES strings"""
        sdf = Path(sims[0].in_path) / f"{sims[0].name}.sdf"
        writer = Chem.SDWriter(str(sdf))

        for sim in sims:
            mol = Chem.MolFromSmiles(sim.smi)
            mol.SetProp("name", sim.name)
            writer.write(mol)
            sim.metadata.prepared_ligand = sdf

        return True

    @staticmethod
    def batch_run(sims: Sequence[Simulation]) -> Optional[list[float]]:
        p_ligand = Path(sims[0].metadata.prepared_ligand)
        ligand_name = p_ligand.stem

        name = f"{Path(sims[0].receptor).stem}_{ligand_name}"

        argv, _, log = SminaRunner.build_argv(
            sims[0].metadata.prepared_ligand,
            sims[0].metadata.prepared_receptor,
            sims[0].center,
            sims[0].size,
            sims[0].ncpu,
            sims[0].metadata.exhaustiveness,
            sims[0].metadata.num_modes,
            sims[0].metadata.energy_range,
            name,
            Path(sims[0].out_path),
            sims[0].metadata.extra,
        )

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            warnings.warn(f'Message: {ret.stderr.decode("utf-8")}', SimulationFailureWarning)

        Ss = SminaRunner.batch_parse_logfile(log)
        if Ss is None:
            s = None
        else:
            s = [utils.reduce_scores(S, sims[0].reduction, k=sims[0].k) for S in Ss]

        for sim, score in zip(sims, s):
            sim.result = Result(
                sim.smi, name, re.sub("[:,.]", "", ray.state.current_node_id()), score
            )

        return s

    @staticmethod
    def batch_parse_logfile(logfile: PathLike) -> Optional[list[np.ndarray]]:
        """parse a Smina-type log file for the scores of the binding modes

        Parameters
        ----------
        logfile : PathLike
            the path to a Vina-type log file

        Returns
        -------
        Optional[list[np.ndarray]]]
            the scores of the docked binding modes in the ordering of the log file. None if no
            scores were parsed or the log file was unparsable
        """
        SCORE_LINE_PATTERN = r"\d\s+\-?\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+"
        pattern = re.compile(SCORE_LINE_PATTERN)

        score_liness = []
        try:
            with open(logfile) as fid:
                for line in fid:
                    while pattern.match(line) is None:
                        line = next(fid)
                    scores_lines = [line]
                    scores_lines.extend(
                        list(takewhile(lambda line: pattern.match(line) is not None, fid))
                    )
                    score_liness.append(scores_lines)
        except StopIteration:
            pass
        except OSError:
            return None

        scoress = [
            np.array([parse_optional_float(line.split()[1]) for line in score_lines])
            for score_lines in score_liness
        ]

        return scoress or None
    
@ray.remote
def batch_prepare_and_run_(sims):
    return SminaRunner.batch_prepare_and_run(sims)

SminaRunner.batch_prepare_and_run_ = batch_prepare_and_run_