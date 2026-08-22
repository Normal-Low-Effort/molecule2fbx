"""PubChem download and RDKit SDF parsing."""

from __future__ import annotations

from typing import Dict, Optional

from .errors import APIError, CIDNotFoundError, InvalidCIDError, No3DConformerError, RDKitError
from .model import Atom, Bond, MoleculeModel, normalize_bond_order


PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def validate_cid(value: object) -> int:
    """Validate a CID supplied by a CLI or another caller."""

    try:
        cid = int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise InvalidCIDError("CID must be a positive integer") from exc
    if cid <= 0 or str(value).strip() != str(cid):
        raise InvalidCIDError("CID must be a positive integer")
    return cid


def _requests_module():
    try:
        import requests
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise APIError("The requests package is required; install the project dependencies first") from exc
    return requests


def _request(url: str, timeout: float):
    requests = _requests_module()
    try:
        return requests.get(
            url,
            timeout=timeout,
            headers={"User-Agent": "molecule2fbx/0.2 (+https://pubchem.ncbi.nlm.nih.gov/)"},
        )
    except requests.RequestException as exc:
        raise APIError(f"Failed to download molecule data: {exc}") from exc


def _cid_exists(cid: int, timeout: float) -> Optional[bool]:
    """Disambiguate a 3D endpoint 404 between an unknown CID and no conformer."""

    url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/Title/JSON"
    try:
        response = _request(url, timeout)
    except APIError:
        return None
    if response.status_code == 200:
        return True
    if response.status_code == 404:
        return False
    return None


def fetch_3d_sdf(cid: int, timeout: float = 30.0) -> str:
    """Download a PubChem 3D SDF, raising a distinct error for each failure mode."""

    cid = validate_cid(cid)
    url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/SDF?record_type=3d"
    response = _request(url, timeout)
    if response.status_code == 200:
        if not response.text.strip():
            raise No3DConformerError("No 3D conformer available")
        return response.text

    body = response.text.lower()
    if response.status_code == 404:
        exists = _cid_exists(cid, timeout)
        if exists is True or any(term in body for term in ("conformer", "3d", "record")):
            raise No3DConformerError("No 3D conformer available")
        raise CIDNotFoundError("CID not found")
    if response.status_code in (400, 422) and any(
        term in body for term in ("conformer", "3d", "record_type")
    ):
        raise No3DConformerError("No 3D conformer available")
    raise APIError(f"Failed to download molecule data (HTTP {response.status_code})")


def _fallback_title(cid: int, timeout: float) -> Optional[str]:
    url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/Title/JSON"
    try:
        response = _request(url, timeout)
    except APIError:
        return None
    if response.status_code != 200:
        return None
    try:
        payload = response.json()
        properties = payload.get("PropertyTable", {}).get("Properties", [])
        title = properties[0].get("Title") if properties else None
    except (AttributeError, IndexError, TypeError, ValueError):
        return None
    return title.strip() if isinstance(title, str) and title.strip() else None


def fetch_compound_properties(cid: int, timeout: float = 30.0) -> Dict[str, object]:
    """Fetch a title and stereochemistry-preserving SMILES for a CID."""

    cid = validate_cid(cid)
    property_sets = ("Title,IsomericSMILES", "Title,SMILES")
    last_status = None
    for property_set in property_sets:
        url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/{property_set}/JSON"
        response = _request(url, timeout)
        last_status = response.status_code
        if response.status_code == 404:
            raise CIDNotFoundError("CID not found")
        if response.status_code != 200:
            continue
        try:
            properties = response.json().get("PropertyTable", {}).get("Properties", [])
            item = properties[0]
        except (AttributeError, IndexError, TypeError, ValueError) as exc:
            raise APIError("PubChem returned invalid compound properties") from exc
        smiles = next(
            (
                item.get(key)
                for key in ("SMILES", "IsomericSMILES", "CanonicalSMILES", "ConnectivitySMILES")
                if isinstance(item.get(key), str) and item.get(key).strip()
            ),
            None,
        )
        return {
            "cid": cid,
            "title": item.get("Title") or f"CID_{cid}",
            "smiles": smiles,
        }
    raise APIError(f"Failed to download molecule properties (HTTP {last_status})")


def parse_sdf(sdf_text: str, cid: int, timeout: float = 30.0) -> MoleculeModel:
    """Parse PubChem's SDF into the Blender-neutral model representation."""

    try:
        from rdkit import Chem
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise RDKitError("The RDKit package is required; install the project dependencies first") from exc

    try:
        mol = Chem.MolFromMolBlock(sdf_text, sanitize=True, removeHs=False, strictParsing=True)
    except Exception as exc:  # RDKit exposes several exception types across versions
        raise RDKitError(f"Could not parse downloaded SDF: {exc}") from exc
    if mol is None:
        raise RDKitError("Could not parse downloaded SDF")
    if mol.GetNumConformers() == 0:
        raise No3DConformerError("No 3D conformer available")

    conformer = mol.GetConformer()
    atoms = []
    for atom in mol.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        atoms.append(
            Atom(
                index=atom.GetIdx(),
                element=atom.GetSymbol(),
                x=float(position.x),
                y=float(position.y),
                z=float(position.z),
            )
        )

    bonds = []
    for bond in mol.GetBonds():
        bonds.append(
            Bond(
                begin=bond.GetBeginAtomIdx(),
                end=bond.GetEndAtomIdx(),
                order=normalize_bond_order(
                    float(bond.GetBondTypeAsDouble()), bond.GetIsAromatic()
                ),
            )
        )

    validated_cid = validate_cid(cid)
    name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
    generic_names = {
        str(validated_cid).casefold(),
        f"cid_{validated_cid}".casefold(),
        f"pubchem cid {validated_cid}".casefold(),
    }
    if not name or name.casefold() in generic_names:
        name = _fallback_title(validate_cid(cid), timeout) or f"CID_{cid}"
    model = MoleculeModel(
        cid=validated_cid,
        name=name,
        atoms=tuple(atoms),
        bonds=tuple(bonds),
        metadata={
            "structure_origin": "pubchem_3d",
            "structure_source": "PubChem PUG REST record_type=3d",
            "formal_charge": int(Chem.GetFormalCharge(mol)),
            "radical_electrons": int(sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())),
            "structure_claim": "PubChem-provided 3D conformer; not labeled as experimental",
        },
    )
    return model


def download_and_parse(cid: int, timeout: float = 30.0) -> MoleculeModel:
    cid = validate_cid(cid)
    return parse_sdf(fetch_3d_sdf(cid, timeout=timeout), cid=cid, timeout=timeout)
