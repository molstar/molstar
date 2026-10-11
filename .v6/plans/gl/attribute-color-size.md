# Plan: Color and size from vertex attributes

## Status

- Not started, to be considered. Shader branches for `dColorType_attribute` / `dSizeType_attribute` (`aColor`, `aSize`) were removed: nothing set these types (`ColorType`/`SizeType` in `geo/geometry/color-data.ts`/`size-data.ts` don't include `'attribute'`), the schema entries were commented-out TODOs, and the branches were never compiled.
- Rule: see [minimum-capabilities.md](minimum-capabilities.md); new attributes count against the 16 vertex attributes (merged cylinders already use all 16) and 8 WebGPU vertex buffers ([attribute-bindings.md](attribute-bindings.md)).

## Idea

Per-vertex color (`vec3`) and size (`float`) supplied directly as attributes instead of texture lookups (`tColor`/`tSize` by instance/group/vertex index).

- Pros: no data texture and index math for geometries that already have per-vertex data (e.g. imported meshes with vertex colors, point clouds); exact values (no 8-bit RGB size packing).
- Cons: more vertex data per instance-independent vertex, no per-instance variation, doesn't combine with group-based themes/overpaint the same way, more attribute locations/vertex buffers (see limits above), and themes produce colors per location (group/instance), not per vertex.

## If implemented

1. Schema: `aColor: AttributeSpec('float32', 3, 0, ['dColorType_attribute'])`, `aSize: AttributeSpec('float32', 1, 0, ['dSizeType_attribute'])`; `'attribute'` back in the `dColorType`/`dSizeType` options and in `ColorType`/`SizeType`.
2. Values: geometries provide the arrays (empty when unused, like `aVertex`); the render item skips buffers whose `usedIf` is known to be false — needs a define-aware variant, since `dColorType` can change at runtime (today only context defines are known when buffers are created).
3. Shaders: restore the branches (`color_vert_params`, `assign_color_varying`, `assign_size`, `dColorType_varying` in `common`); declarations are generated.
4. Themes: a color/size theme kind that provides per-vertex values, or a geometry-level option.
5. Tests: limits test (attribute locations/vertex buffers with the new attributes), visual cases.
