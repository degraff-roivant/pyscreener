from dataclasses import asdict
from typing import Dict, Optional

from colorama import init, Fore, Style

from ..exceptions import InvalidEnvironmentError, UnsupportedSoftwareError
from .sim import Simulation
from .metadata import SimulationMetadata
from .result import Result
from .runners import DockingRunner
from .screen import DockingVirtualScreen
from .utils import ScreenType
from . import vina, smina, dock

init(True)

VALID_SOFTWARE = {"vina", "qvina", "psovina", "smina", "dock", "dock6", "ucsfdock"}

def build_metadata(software: str, metadata: Optional[Dict] = None) -> SimulationMetadata:
    software_ = software.lower()
    metadata = metadata or {}

    if software_ in ("vina", "qvina", "psovina"):
        d_md = asdict(vina.VinaMetadata())
        d_md.update((k, metadata[k]) for k in d_md.keys() & metadata.keys())

        return vina.VinaMetadata(**d_md)

    if software_ == "smina":
        d_md = asdict(smina.SminaMetadata())
        d_md.update((k, metadata[k]) for k in d_md.keys() & metadata.keys())

        return smina.SminaMetadata(**d_md)

    if software_ in ("dock", "dock6", "ucsfdock"):
        d_md = asdict(dock.DOCKMetadata())
        d_md.update((k, metadata[k]) for k in d_md.keys() & metadata.keys())

        return dock.DOCKMetadata(**d_md)

    raise UnsupportedSoftwareError(
        f'Unrecognized screen type! got: "{software}. Expected one of {VALID_SOFTWARE}."'
    )


def get_runner(software: str) -> DockingRunner:
    software_ = software.lower()

    if software_ == "smina":
        return smina.SminaRunner

    if software_ in ("vina", "qvina", "psovina"):
        return vina.VinaRunner

    if software_ in ("dock", "dock6", "ucsfdock"):
        return dock.DOCKRunner

    raise UnsupportedSoftwareError(
        f'Unrecognized screen type! got: "{software}. Expected one of {VALID_SOFTWARE}."'
    )


def check_env(software, metadata: Optional[Dict] = None):
    print(f'Checking environment and metadata for "{software}" screen')
    try:
        valid_env = False
        print("  Checking PATH and environment variables ...", end=" ")
        runner = get_runner(software)
        print(Style.BRIGHT + Fore.GREEN + "PASS")
        valid_env = True
        print("  Validating metadata ... ", end=" ")
        metadata = build_metadata(software, metadata)
        runner.check_environment(metadata)
        print(Style.BRIGHT + Fore.GREEN + "PASS")
    except (InvalidEnvironmentError, UnsupportedSoftwareError):
        print(Style.BRIGHT + Fore.RED + "FAIL")
        if not valid_env:
            print("Environment not set up properly!", end=" ")
        else:
            print("Invalid metadata supplied!", end=" ")
        print("See the exception for more details", flush=True)
        raise

    print("Environment is properly set up!")


def virtual_screen(software: str, *args, **kwargs) -> DockingVirtualScreen:
    return DockingVirtualScreen(get_runner(software), *args, **kwargs)
