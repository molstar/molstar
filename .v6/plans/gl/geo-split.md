# Plan: the non-GL part of geo

## Status

- **Current**: step 1 done (the marching cubes computation in core); step 2 deferred to step 4's decision (where the
  geometry data goes decides where the primitives go: one move, not two); step 3 proposed, step 4 needs its design
  decision first.
- Goal: geometry math usable without graphics (core, io, model, servers, the extensions' algorithms), as
  `core/math/geometry` already holds the molecular surface, the Gaussian density on the CPU and the boundaries. Graphics
  keeps what renders (values, schemas, themes) and what is a rendering concept though GL-free.

## Inventory

`geo`'s files (outside `_test`) by what they reach, directly or through other `geo` files:

| Part | Files | Reaches | Destination |
| --- | --- | --- | --- |
| Primitives (`geo/primitive/*`: polyhedra, prism, torus, ... as vertices, normals, indices, `PrimitiveBuilder`) | 20 | core | core (step 2) |
| Array math (`geo/util.ts`: array transforms, vertex normals, radii), Laplacian smoothing | 2 | core | core, where array-based (step 2) |
| Marching cubes computation, tables, estimate | 3 | core | done (step 1): `core/math/geometry/marching-cubes` |
| Rendering concepts, GL-free (animation, interior, picking ids, transform data) | 4 | core | stay: values and parameters of rendering |
| Location iterators, streamlines | 2 | model | stay: themes' locations, model properties |
| Geometry types and their builders (`Mesh`, `Lines`, `Spheres`, ... with `MeshBuilder`, `LinesBuilder`, the box, cylinder, tube, ribbon, sheet builders, the marching cubes builders) | ~30 | gl, theme | step 4: the data split from the rendering, then decided |
| Value data (color, size, marker, overpaint, ...), color smoothing, GPU density, texture mesh, shapes | ~15 | gl, theme, repr | stay: rendering |

## Steps

1. **Marching cubes** (done): `computeMarchingCubes(ctx, builder, params)`, `MarchingCubesBuilder`, the tables and the
   estimate in `core/math/geometry/marching-cubes`; graphics keeps the mesh and lines builders and
   `computeMarchingCubesMesh`/`computeMarchingCubesLines`. The GPU marching cubes import the tables from core (S7's
   `gl -> geo` edge gone).
2. **Primitives and array math to core**: `geo/primitive/*` and the array functions of `geo/util.ts` (and the Laplacian
   smoothing where it works on arrays, not on a `Mesh`) into `core/math/geometry` (a folder of their own: `primitives`
   is taken there by `Box3D`, `Sphere3D` and the like, so e.g. `core/math/geometry/polyhedra`). Graphics imports them;
   its mesh builders keep adding primitives to meshes.
3. **A guard**: a source test that the moved modules and the rest of `core` import nothing of graphics (with
   `core-gpu-compute.md`'s step 1, one test).
4. **Geometry data from rendering**: each geometry module (`mesh.ts`, `lines.ts`, ...) holds the data (buffers in
   `ValueCell`s, counts, bounding sphere, group mapping) and, in `Mesh.Utils` and the like, its rendering (params,
   `createValues`, `updateValues`, the renderable state, from `gl` schemas and a theme). Split them per geometry: a data
   module importing core only, a rendering module with the utils. Then the builders reach core only, and the decision
   below can be made.

## Decisions open

- **Where the geometry data and builders go** after step 4: into core (`core/math/geometry/mesh` and the like, the
  builders with them: any package can build a mesh), a package of their own between core and graphics
  (`@molstar/geometry`: the geometry data, builders and primitives, maybe the surface algorithms; core keeps no
  graphics data format, graphics renders it), or a GL-free folder of graphics (`geo/data`, nothing moves out of the
  package, the boundary held by a test). The cost of core: a `Mesh` is render-ready geometry (its `ValueCell` buffers
  are what the renderables bind), so core would hold a graphics data format. A package for the primitives alone buys
  nothing: everything using them depends on core already.
- **The primitives' place** (step 2): with the geometry data where that goes; in core, a folder of their own in
  `core/math/geometry` (`polyhedra` or `shapes`; not `core/geometry`, which would split geometry in two).

## Why not more

The value data, color smoothing, GPU density and texture mesh render: they stay. The rendering concepts that are GL-free
today (animation, interior, picking ids, transform data) are values and parameters of renderables: in core they would
be graphics' vocabulary without graphics.
