## Graphics: faster mesh building

- Mesh primitives (dashed interaction cylinders, mesh spheres, carbohydrate symbols, `addMesh`) and `Mesh.transform`
  read and transform their arrays inline instead of through `Vec3.fromArray` and the `Vec3` transforms, which are
  megamorphic app-wide; identical output. [Breaking] `Primitive` arrays are `Float32Array` / `Uint32Array`, see the
  [API change ledger](../plans/breaking-v6-changes.md)
