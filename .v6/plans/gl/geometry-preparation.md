# Plan: Geometry preparation before values

## Status

- **Current**: done (steps 1–4 except geo-export). `GeometryRequest` (`geo/geometry/request.ts`, `fromTheme` from
  `isPositionColorType`/`isPositionSizeType` next to `createColors`/`createSizes`), `GeometryUtils.prepare`, the texture
  mesh's `prepare` (reads after `meta.extracted`, which an indirectly drawn mesh sets to its count check and any
  re-extraction), `Visual.withPreparedGeometry` in the units, complex, volume and particle visuals;
  `TextureMesh.withVertexData`, `isVertexDataNeeded` and the exported `readVertexData` are gone. Unit tests
  (`geo/geometry/_test/request.test.ts`, `repr/_test/prepared-geometry.test.ts`); the visual suite passes on `webgl2`,
  `webgl1`, `webgpu` and `webgpu-core` (`gpu-paths/vertex-colors-gpu-surface` on WebGL2 and WebGPU). The shape
  representation stays as it is: its geometry comes from `getShape` on the CPU and it has no graphics context to prepare
  with. geo-export keeps its own readback (it needs groups and indices as well, a request for them when it moves).

## Where things stand

- **What it solves**: a GPU surface (`TextureMesh`) keeps its vertices in textures. The position iterator, which color
  and size themes use to sample at vertex positions (`granularity` `vertex`/`vertexInstance`, e.g. coloring a surface by
  an external volume; `createColors` and `createSizes` are its only consumers), needs them on the CPU.
  `TextureMesh.withVertexData` reads them back asynchronously (`readVertexData` into `meta.vertexData`) before a visual
  creates or updates its values, and only when the theme needs them (`isVertexDataNeeded`, which checks the color theme
  only: a vertex-granularity size theme gets no positions; harmless today as texture meshes have no size values, not for
  future kinds with GPU positions).
- **How it is wired**: the units, complex and volume visual builders call `TextureMesh.withVertexData` directly, each in
  three branches (geometry created synchronously, created asynchronously, kept) that wrap `update` and `finalize`. The
  particle visual and the shape representation don't call it, so a texture mesh there gets `NullLocation` positions.
- **Invalidation is spread out**: `TextureMesh.create` and `setVertexCount` clear `meta.vertexData`; re-extraction
  handlers (`meta.reextracted`, after an indirectly drawn mesh's vertex count was read back) re-run color smoothing but
  not the read back.
- **Similar needs elsewhere**: geo-export reads texture meshes back on its own (`readTextureData`); planned GPU paths
  (positions of spheres or cylinders from compute, [core-gpu-compute.md](core-gpu-compute.md)) would need the same
  CPU-side data for iterators, bounding spheres or export.

## Design

### A preparation step per geometry kind

`GeometryUtils` gets an optional step that readies geometry-derived data the values need:

```ts
interface GeometryUtils<G> {
  /**
   * Readies data of `geometry` that creating or updating values for `theme` needs (e.g. GPU data read back for the
   * position iterator); undefined when nothing is needed (synchronous kinds omit it)
   */
  prepare?(geometry: G, request: GeometryRequest, graphics: GraphicsContext | undefined): Promise<void> | undefined;
}

/** plain flags, `prepare` sees no themes */
interface GeometryRequest {
  /** the position iterator is used (vertex-granularity color or size themes) */
  positions: boolean;
}
```

- **Request in the geometry layer**: `GeometryRequest` and its derivation `GeometryRequest.fromTheme(theme)` live in the
  geometry layer, next to `createColors` (`color-data.ts`) and `createSizes` (`size-data.ts`), which decide which
  granularities use the position iterator: `positions` is true for a `vertex`/`vertexInstance` color or size theme. The
  geometry layer takes themes already (`GeometryUtils.createValues`, `createColors`, `createSizes`), so this adds no
  dependency; the visuals pass the theme through and don't know which themes need what. Overpaint, transparency,
  emissive and substance work per group (their smoothing runs on the GPU), clipping and marking per group: no positions,
  no flags. A further need (normals and indices for export) is a further flag, without new entry points.
- The texture mesh implements `prepare` with today's `isVertexDataNeeded`/`readVertexData`; other kinds omit it.
- Returning undefined when nothing is to do keeps synchronous updates synchronous (the synchronous runtime of headless
  rendering and tests).

### One update path in the visuals

A shared helper replaces the three branches in each builder:

```ts
/** `geometry` (created, possibly asynchronously, or kept) prepared for `theme`, then `update` */
function withPreparedGeometry<G>(
  geometry: G | Promise<G>,
  utils: GeometryUtils<G>,
  theme: Theme,
  graphics: GraphicsContext | undefined,
  update: (geometry: G) => void,
): Promise<void> | void;
```

- Units, complex, volume and particle visuals use it (`Visual.withPreparedGeometry`), so every geometry kind gets
  prepared wherever a visual draws it, including the particle visual that missed it.
- `TextureMesh.withVertexData` and `isVertexDataNeeded` go away; `readVertexData` stays as the texture mesh's
  implementation detail.

### Invalidation owned by the geometry kind

- The kind clears its prepared data when the geometry changes (texture mesh: `create`, `setVertexCount`) and `prepare`
  reads again on the next update. An indirectly drawn mesh may be extracted again once its count is read back: `prepare`
  waits for that (`meta.extracted`), so the values are created from the final vertices.
- A backend switch needs nothing extra: `prepare` runs with the visual's current context.

## Steps

1. `GeometryRequest` with `fromTheme` (color and size), `GeometryUtils.prepare`, the texture mesh's implementation; unit
   tests with a geometry kind whose `prepare` resolves later (values created after it, synchronous kinds stay
   synchronous) and of `fromTheme` for vertex-granularity size themes.
2. `withPreparedGeometry`, used by the units, complex and volume visuals in place of `TextureMesh.withVertexData`
   (behavior unchanged: `gpu-paths/vertex-colors-gpu-surface` passes on all profiles).
3. The particle visual and the shape representation on the same path; a case for a texture mesh with vertex colors in
   one of them if a representation produces one.
4. Invalidation in the texture mesh (re-extraction requests the read back again); geo-export reads through `prepare` (a
   request for positions, normals and indices) instead of its own readback, if that fits.
