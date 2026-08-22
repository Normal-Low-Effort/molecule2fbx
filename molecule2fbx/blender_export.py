"""Invoke Blender's Python API in a separate, headless process."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from .errors import BlenderExportError, BlenderNotFoundError
from .model import MoleculeModel, validate_model


def find_blender(explicit: Optional[str] = None) -> str:
    """Find Blender from an explicit path, environment variable, or PATH."""

    candidates = []
    if explicit:
        candidates.append(explicit)
    env_value = os.environ.get("BLENDER_EXECUTABLE")
    if env_value:
        candidates.append(env_value)
    candidates.extend(["blender", "blender.exe"])

    for candidate in candidates:
        path = Path(candidate)
        if path.is_file():
            return str(path.resolve())
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    raise BlenderNotFoundError(
        "Blender executable not found. Install Blender or pass --blender <path>."
    )


def _worker_path() -> Path:
    return Path(__file__).with_name("blender_worker.py")


def _error_tail(result: subprocess.CompletedProcess[str]) -> str:
    output = (result.stderr or result.stdout or "").strip()
    if not output:
        return "no Blender output"
    lines = output.splitlines()
    return "\n".join(lines[-12:])


def export_model(
    model: MoleculeModel,
    output_path: Path,
    blender_executable: str,
    timeout: float = 180.0,
) -> Path:
    """Serialize model data, run Blender, and verify that the FBX was created."""

    validate_model(model)
    worker = _worker_path()
    if not worker.is_file():
        raise BlenderExportError(f"Bundled Blender worker not found: {worker}")

    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="molecule2fbx-") as temp_dir:
        data_path = Path(temp_dir) / "molecule.json"
        data_path.write_text(
            json.dumps(model.to_dict(), ensure_ascii=False, separators=(",", ":")),
            encoding="utf-8",
        )
        command = [
            blender_executable,
            "--background",
            "--factory-startup",
            "--python",
            str(worker),
            "--",
            "--data-file",
            str(data_path),
            "--output",
            str(output_path),
        ]
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=timeout,
                check=False,
            )
        except FileNotFoundError as exc:
            raise BlenderNotFoundError(
                f"Blender executable could not be started: {blender_executable}"
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise BlenderExportError(f"Blender export timed out after {timeout:g} seconds") from exc
        except OSError as exc:
            raise BlenderExportError(f"Could not start Blender: {exc}") from exc

    if result.returncode != 0:
        raise BlenderExportError(
            f"Blender failed to export FBX (exit code {result.returncode})\n{_error_tail(result)}"
        )
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise BlenderExportError("Blender finished without creating an FBX file")
    return output_path
