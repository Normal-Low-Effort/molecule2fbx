"""Electronic-state validation and conservative metal detection."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple

from .errors import ConfigurationError, UnsupportedElementError


_SYMBOLS = (
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
    "Ts", "Og",
)
_ATOMIC_NUMBERS = {symbol: number for number, symbol in enumerate(_SYMBOLS) if symbol}
_METALS = {
    3, 4, 11, 12, 13, 19, 20,
    *range(21, 32),
    37, 38,
    *range(39, 51),
    55, 56,
    *range(57, 85),
    87, 88,
    *range(89, 113),
}


@dataclass(frozen=True)
class ElectronicState:
    electrons: int
    charge: int
    multiplicity: int
    metal_elements: Tuple[str, ...]


def validate_electronic_state(
    elements: Iterable[str],
    charge: int,
    multiplicity: int,
) -> ElectronicState:
    atomic_numbers = []
    for element in elements:
        number = _ATOMIC_NUMBERS.get(element)
        if number is None:
            raise UnsupportedElementError(f"Unknown or unsupported element symbol: {element}")
        atomic_numbers.append(number)
    electrons = sum(atomic_numbers) - charge
    if electrons <= 0:
        raise ConfigurationError(
            f"Charge {charge} leaves a non-positive electron count ({electrons})"
        )
    if multiplicity < 1:
        raise ConfigurationError("Multiplicity must be a positive integer")
    unpaired = multiplicity - 1
    if unpaired > electrons or (electrons - unpaired) % 2:
        raise ConfigurationError(
            "Charge/multiplicity is inconsistent with the electron count: "
            f"electrons={electrons}, charge={charge}, multiplicity={multiplicity}"
        )
    metals = tuple(sorted({_SYMBOLS[z] for z in atomic_numbers if z in _METALS}))
    return ElectronicState(electrons, charge, multiplicity, metals)


def validate_metal_policy(
    state: ElectronicState,
    *,
    method: str,
    allow_metals: bool,
    basis_explicit: bool,
    functional_explicit: bool,
    charge_explicit: bool,
    multiplicity_explicit: bool,
) -> None:
    if not state.metal_elements:
        return
    element_list = ", ".join(state.metal_elements)
    if not allow_metals:
        raise UnsupportedElementError(
            f"Metal-containing molecule detected ({element_list}). The organic-molecule defaults "
            "will not be applied automatically; review the basis/ECP and spin state, then pass "
            "--allow-metals with explicit calculation settings."
        )
    required = [
        (basis_explicit, "--basis"),
        (charge_explicit, "--charge"),
        (multiplicity_explicit, "--multiplicity"),
    ]
    if method == "dft":
        required.append((functional_explicit, "--functional"))
    missing = [option for supplied, option in required if not supplied]
    if missing:
        raise ConfigurationError(
            "Metal calculations require explicit settings: " + ", ".join(missing)
        )
