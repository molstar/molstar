## Graphics: faster mesh building

- [Breaking] `Primitive` arrays are `Float32Array` / `Uint32Array` (all built-in producers already made them), see the
  [API change ledger](../plans/breaking-v6-changes.md)
- Mesh primitives (dashed interaction cylinders, mesh spheres, carbohydrate symbols, `MeshBuilder.addMesh`), the
  geometry transforms (`Mesh.transform`, lines, cylinders, points), `Mesh.computeNormals` and the cartoon curve and
  builders read their arrays inline instead of through `Vec3.fromArray` and the `Vec3` transforms, which go megamorphic
  app-wide; output is bit-identical. Node, 300k capped cylinders (30M vertices): `addPrimitive` 469 → 305 ms, 966 →
  312 ms once `Vec3` has seen several array kinds (as in the app); transforming 15M positions and normals 83 → 46 ms
  (2.6 s → 47 ms polymorphic); vertex normals for 6M vertices 47 → 29 ms (149 → 29 ms polymorphic) (#1990)
- The cylinder primitive has no UV duplicates: no seam column in a full turn, one center vertex per cap (16 capped
  segments: 66 vertices instead of 100), same triangles. Interactions at 186k atoms: 998 → 373 ms, 623 → 452 MB of
  buffers (#1990)
