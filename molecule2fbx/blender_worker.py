"""Blender-only worker script used by :mod:`molecule2fbx.blender_export`.

This file is executed by Blender and intentionally imports ``bpy`` only in that
process. It is kept as a normal package file so an installed CLI can locate it.
"""

from __future__ import print_function

import json
import os
import sys
import traceback

import bpy
from mathutils import Vector


ATOM_COLORS = {
    "H": (0.95, 0.95, 0.95),
    "C": (0.12, 0.12, 0.12),
    "N": (0.08, 0.20, 0.85),
    "O": (0.85, 0.05, 0.05),
    "F": (0.10, 0.45, 0.10),
    "Cl": (0.55, 0.80, 0.10),
    "S": (0.95, 0.80, 0.05),
    "P": (0.95, 0.45, 0.08),
    "Br": (0.55, 0.20, 0.08),
    "I": (0.35, 0.08, 0.55),
}

ATOM_RADII = {
    "H": 0.25,
    "C": 0.38,
    "N": 0.33,
    "O": 0.30,
    "F": 0.29,
    "Cl": 0.36,
    "S": 0.37,
    "P": 0.39,
    "Br": 0.38,
    "I": 0.40,
}


def _arguments():
    try:
        separator = sys.argv.index("--")
        args = sys.argv[separator + 1 :]
    except ValueError:
        args = []
    values = {}
    index = 0
    while index < len(args):
        key = args[index]
        if key in ("--data-file", "--output") and index + 1 < len(args):
            values[key[2:]] = args[index + 1]
            index += 2
        else:
            index += 1
    if not values.get("data-file") or not values.get("output"):
        raise ValueError("Blender worker requires --data-file and --output")
    return values


def _material(name, color):
    material = bpy.data.materials.get(name) or bpy.data.materials.new(name=name)
    material.use_nodes = True
    material.diffuse_color = (*color, 1.0)
    nodes = material.node_tree.nodes
    shader = nodes.get("Principled BSDF")
    if shader:
        shader.inputs["Base Color"].default_value = (*color, 1.0)
        if "Roughness" in shader.inputs:
            shader.inputs["Roughness"].default_value = 0.32
    return material


def _add_uv_sphere(atom, materials, root):
    element = atom["element"]
    location = Vector(atom["position"])
    radius = ATOM_RADII.get(element, 0.32)
    bpy.ops.mesh.primitive_uv_sphere_add(
        segments=32,
        ring_count=20,
        radius=radius,
        location=location,
    )
    obj = bpy.context.object
    obj.name = "Atom_{0}_{1}".format(atom["index"], element)
    obj.data.materials.append(materials.get(element, materials["__default__"]))
    obj["element"] = element
    obj["atom_index"] = atom["index"]
    obj.parent = root
    for polygon in obj.data.polygons:
        polygon.use_smooth = True
    return obj


def _cylinder_between(start, end, radius, material, name, root):
    vector = end - start
    length = vector.length
    if length < 1.0e-7:
        return None
    midpoint = (start + end) * 0.5
    bpy.ops.mesh.primitive_cylinder_add(
        vertices=24,
        radius=radius,
        depth=length,
        location=midpoint,
    )
    obj = bpy.context.object
    obj.name = name
    obj.rotation_mode = "QUATERNION"
    obj.rotation_quaternion = vector.to_track_quat("Z", "Y")
    obj.data.materials.append(material)
    obj.parent = root
    for polygon in obj.data.polygons:
        polygon.use_smooth = True
    return obj


def _offsets(start, end, order):
    if order <= 1:
        return [Vector((0.0, 0.0, 0.0))]
    axis = (end - start).normalized()
    reference = axis.cross(Vector((0.0, 0.0, 1.0)))
    if reference.length < 1.0e-7:
        reference = axis.cross(Vector((0.0, 1.0, 0.0)))
    reference.normalize()
    spacing = 0.18
    if order == 2:
        return [-reference * (spacing * 0.5), reference * (spacing * 0.5)]
    return [-reference * spacing, Vector((0.0, 0.0, 0.0)), reference * spacing]


def _add_bond(bond, atoms_by_index, materials, root):
    begin_atom = atoms_by_index[bond["begin"]]
    end_atom = atoms_by_index[bond["end"]]
    start = Vector(begin_atom["position"])
    end = Vector(end_atom["position"])
    order = int(bond["order"])
    radius = 0.11 if order == 1 else 0.085
    begin_material = materials.get(begin_atom["element"], materials["__default__"])
    end_material = materials.get(end_atom["element"], materials["__default__"])
    for offset_index, offset in enumerate(_offsets(start, end, order)):
        shifted_start = start + offset
        shifted_end = end + offset
        midpoint = (shifted_start + shifted_end) * 0.5
        prefix = "Bond_{0}_{1}_order{2}_{3}".format(
            bond["begin"], bond["end"], order, offset_index
        )
        _cylinder_between(
            shifted_start, midpoint, radius, begin_material, prefix + "_A", root
        )
        _cylinder_between(midpoint, shifted_end, radius, end_material, prefix + "_B", root)


def build_scene(data):
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.scale_length = 1.0

    materials = {
        element: _material("Element_" + element, color)
        for element, color in ATOM_COLORS.items()
    }
    materials["__default__"] = _material("Element_Default", (0.55, 0.55, 0.55))

    root = bpy.data.objects.new("Molecule_Metadata", None)
    scene.collection.objects.link(root)
    root["molecule_name"] = data["name"]
    if data.get("cid") is not None:
        root["pubchem_cid"] = int(data["cid"])
    for key, value in data.get("metadata", {}).items():
        if value is None:
            continue
        if isinstance(value, (str, int, float, bool)):
            root[key] = value
        else:
            root[key] = json.dumps(value, ensure_ascii=False, sort_keys=True)

    atoms = data["atoms"]
    atoms_by_index = {atom["index"]: atom for atom in atoms}
    for atom in atoms:
        _add_uv_sphere(atom, materials, root)
    for bond in data["bonds"]:
        _add_bond(bond, atoms_by_index, materials, root)

    scene["molecule_name"] = data["name"]
    if data.get("cid") is not None:
        scene["pubchem_cid"] = int(data["cid"])
    scene["structure_origin"] = data.get("metadata", {}).get("structure_origin", "unknown")


def export_fbx(output_path):
    output_path = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    try:
        bpy.ops.export_scene.fbx(
            filepath=output_path,
            use_selection=False,
            object_types={"MESH", "EMPTY"},
            path_mode="AUTO",
            bake_anim=False,
            add_leaf_bones=False,
            use_custom_props=True,
        )
    except TypeError:
        # Keep compatibility with older Blender FBX exporters.
        bpy.ops.export_scene.fbx(
            filepath=output_path,
            use_selection=False,
            object_types={"MESH", "EMPTY"},
            path_mode="AUTO",
            bake_anim=False,
            add_leaf_bones=False,
        )
    if not os.path.isfile(output_path) or os.path.getsize(output_path) == 0:
        raise RuntimeError("FBX exporter did not create an output file")


def main():
    arguments = _arguments()
    with open(arguments["data-file"], "r", encoding="utf-8") as handle:
        data = json.load(handle)
    build_scene(data)
    export_fbx(arguments["output"])
    print("MOLECULE2FBX_OUTPUT=" + os.path.abspath(arguments["output"]))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("molecule2fbx Blender worker error: {0}".format(exc), file=sys.stderr)
        traceback.print_exc()
        raise
