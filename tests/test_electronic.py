import pytest

from molecule2fbx.electronic import validate_electronic_state, validate_metal_policy
from molecule2fbx.errors import ConfigurationError, UnsupportedElementError


def test_validates_electron_parity_before_quantum_execution():
    state = validate_electronic_state(["O", "H", "H"], charge=0, multiplicity=1)
    assert state.electrons == 10
    with pytest.raises(ConfigurationError, match="inconsistent"):
        validate_electronic_state(["O", "H", "H"], charge=0, multiplicity=2)


def test_metal_defaults_are_blocked():
    state = validate_electronic_state(["Fe", "C", "H", "H"], charge=0, multiplicity=1)
    with pytest.raises(UnsupportedElementError, match="Metal-containing"):
        validate_metal_policy(
            state,
            method="dft",
            allow_metals=False,
            basis_explicit=False,
            functional_explicit=False,
            charge_explicit=False,
            multiplicity_explicit=False,
        )


def test_metal_override_requires_explicit_settings():
    state = validate_electronic_state(["Fe"], charge=0, multiplicity=1)
    with pytest.raises(ConfigurationError, match="--basis"):
        validate_metal_policy(
            state,
            method="dft",
            allow_metals=True,
            basis_explicit=False,
            functional_explicit=True,
            charge_explicit=True,
            multiplicity_explicit=True,
        )
