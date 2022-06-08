from enum import auto
from typing import Iterable

from pyscreener.utils import AutoName


class ScreenType(AutoName):
    DOCK = auto()
    VINA = auto()


def parse_optional_float(s: str):
    try:
        return float(s)
    except ValueError:
        return None
        