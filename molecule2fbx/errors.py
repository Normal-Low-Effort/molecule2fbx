"""Domain-specific exceptions used by molecule2fbx."""


class Molecule2FBXError(Exception):
    """Base class for expected, user-facing failures."""


class InvalidCIDError(Molecule2FBXError):
    """The command-line CID is not a positive integer."""


class CIDNotFoundError(Molecule2FBXError):
    """PubChem does not know the requested CID."""


class No3DConformerError(Molecule2FBXError):
    """The compound exists but PubChem has no 3D conformer for it."""


class APIError(Molecule2FBXError):
    """PubChem could not be reached or returned an unexpected response."""


class RDKitError(Molecule2FBXError):
    """RDKit is unavailable or could not parse the downloaded structure."""


class InvalidSMILESError(Molecule2FBXError):
    """The supplied SMILES could not be parsed."""


class ConfigurationError(Molecule2FBXError):
    """The requested calculation settings are inconsistent or unsafe."""


class ForceFieldError(Molecule2FBXError):
    """RDKit could not generate or optimize an initial conformer."""


class QuantumBackendNotFoundError(Molecule2FBXError):
    """The requested external quantum chemistry backend is unavailable."""


class QuantumCalculationError(Molecule2FBXError):
    """The quantum chemistry process failed for an unexpected reason."""


class SCFNotConvergedError(QuantumCalculationError):
    """The electronic self-consistent-field calculation did not converge."""


class GeometryNotConvergedError(QuantumCalculationError):
    """Geometry optimization ended without satisfying convergence criteria."""


class UnsupportedElementError(QuantumCalculationError):
    """The selected defaults must not be applied to the molecule's elements."""


class BlenderNotFoundError(Molecule2FBXError):
    """No Blender executable could be found."""


class BlenderExportError(Molecule2FBXError):
    """Blender failed to create a valid FBX file."""
