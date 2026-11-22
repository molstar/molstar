# Plan: Gaussian splats

## Status

- **Current**: steps 1 to 9 done: the splat data and the data layers' page bits by the schema's own layers; the kind
  (`GaussianSplats.Utils`, the schema, `gaussian-splats` programs in GLSL and WGSL), the draw fact and the
  `gaussianSplats` pass in each transparency mode and illumination; the order back to front on the CPU (`tSplatOrder`);
  PLY files of both layouts; picking and marking where at least half opaque, group bounds, layers; the order on the GPU
  (WebGPU); view-dependent color (spherical harmonics); `.splat`, `.spz` and `.sog` files (read by layers: `io`,
  `model`, `graphics`). Checked on WebGPU (Dawn: the projected Gaussian's widths within 1% of the CPU's, two splats in
  order from either side), in Chrome (visual cases `geometry/gaussian-splats/**` and `transparency/*/splats`: WebGL2 and
  WebGPU alike, by eye) and in the Viewer (PlayCanvas's compressed guitar and biker on WebGL2 and WebGPU alike,
  upright; WebGL1 as WebGL2, the same splat picked). Standard-layout files of degree 3:
  Postshot's cactus (139k splats) on WebGL2 and WebGPU alike, the strawberry (1.5 M: loaded in 1.7 s, its view-dependent
  color 2.3 s more) on WebGPU. Next: the open items. A v5 prototype (a stash, `gaussian-splats`, on 5.10.1) is the
  starting point for the math and the PLY decoding, not for the wiring (see Salvage).
- Related: [primitive-expansion.md](primitive-expansion.md) (per-item data, quads from the vertex id),
  [storage-buffers.md](storage-buffers.md) (`DataSpec`, data pages), [modular-passes.md](modular-passes.md) §4 (new
  renderer passes, the transmission pass as the template), [draw-table.md](draw-table.md) (draw facts),
  [pipeline-readiness.md](pipeline-readiness.md) (passes collected), [core-gpu-compute.md](core-gpu-compute.md) (no
  sort on WebGL), [neutral-conventions.md](neutral-conventions.md).

## Question

Render 3D Gaussian splatting scenes (3DGS: anisotropic Gaussians with a color and an opacity, millions of them) as a
geometry kind next to the molecular ones, on both backends: loaded from PLY (the 3DGS layout and the compressed
PlayCanvas/SuperSplat layout), colored and picked as a shape, composited with the rest of the scene.

## Current state

- **Reader**: `@molstar/io`'s PLY parser reads ASCII and binary (both endians) from bytes; binary tables are read a
  value at a time through a `DataView` (fine for meshes, slow for millions of rows of ~62 properties). `shapeFromPly`
  (`graphics/src/formats/shape/ply.ts`) requires a `face` element: a splat PLY fails with 'missing face element'.
  `PlyProvider` (`plugin/core/src/state/formats/shape/ply.ts`) owns the `ply` extension.
- **Geometry**: quad kinds (points, spheres) hold one `tPositionGroup` texel an item and expand 6 vertex ids in the
  vertex shader (`assign_quad_position_group`, shared quad indices, `QuadItems`). No geometry orders its items per
  view.
- **Data layers**: a `DataSpec` reaches WebGPU only if named in `DataLayers` (`gl/webgpu/data-layers.ts`). The page
  define `dDataLayerPages` has 2 bits a layer in the order of the schema's layers (`getSchemaDataLayers`): 16 layers a
  schema (`MaxSchemaDataLayers`, checked), any number in the list (it held 16, in the list's order, before step 1).
  WebGL has no such limit.
- **Transparency**: `DrawPass` runs one of WBOIT, DPOIT or blended per eye (`canvas3d/passes/draw.ts`), illumination a
  copy of the same (`passes/illumination.ts`); passes are rows of `PassDraws` (`gl/renderer/pass-draws.ts`) selected by
  `DrawFact`s (`gl/draw-table.ts`, bits 0-27 used). Objects are drawn in program order, never back-to-front.
- **Per-view orders**: instance expansion sorts copies by distance on the CPU into `tCopyOrder`, kept for up to 4
  views (`CopyOrderViews`, `gl/instance-expansion.ts`) and uploaded lazily before a draw (`syncCopyOrder`). That is the
  pattern for a splat order. `@molstar/core` has `radixSortIndices` (`data/util/sort.ts`).
- **Compute**: WebGPU only (`ComputeContext`), with a prefix sum, atomics and indirect dispatch, no sort; a cull batch
  per eye before the passes (`Renderable.prepareCull`). WebGL's fragment-pass compute cannot sort.

## Proposal

### Data

Per splat, once, in the layout the shaders read:

| Layer            | Format             | Content                                                                                  | Bytes |
| ---------------- | ------------------ | ---------------------------------------------------------------------------------------- | ----- |
| `tPositionGroup` | `rgb32f-i32`       | center, group (shared with points and spheres: bounds, group mapping, position iterator) | 16    |
| `tSplat`         | `rgba16i`, 2 a splat | rotation x, y, z snorm16 and w unorm16; log-scales x, y, z in fixed point, opacity unorm16 | 16    |

- **Scale and rotation, not covariance**: the vertex shader builds `Σ = R S Sᵀ Rᵀ` (a few dozen flops a vertex); 16
  bytes instead of 28 for six f32 plus opacity, and f16 covariance underflows for small splats (scale² of 1e-3 Å).
- **Encoding** (`splat-data.ts`): `rgba16i` as spheres' sites (no new format; WebGL1, without integer textures, has
  float pages, the integers converted where written as `rgb32f-i32`'s groups: `convertInt16ToFloat`): w kept non-negative
  (`q` and `-q` one rotation) so the unsigned fourth holds it; the natural logs of the scales (as 3DGS stores them) in
  steps of 1/1024 (logs within ±32, relative steps of 0.1%: uniform precision at any size, unlike f16 scales).
- **Color** is the theme's per group (`tColor`, 'group' granularity for the shape's own colors): the PLY's DC term is
  the shape's color function, as in the prototype, so overpaint, uniform color and the other themes apply. A splat
  representation's default color theme is the shape's.
- **Groups**: one a splat by default (a splat picked and highlighted alone); a 'none' grouping (one group) for scenes
  that only display, which makes markers one value.
- **The 17th layer**: the page bits by the schema's layers, not by the list's (decided in step 1, over a storage
  buffer for `tSplat` on WebGPU alone): a splat schema has 11 layers, the list any number.

### Order

- **Key**: the view z of a splat's center (`a x + b y + c z` by the z row of model-view, object and instance transform;
  perspective and orthographic alike), quantized to 16 bits over the instance's range, back to front (the most negative
  first) by a counting sort (`countingSortIndices` in `@molstar/core`: stable, one pass; 1 M splats 8-10 ms in Node, 5 M
  42-47 ms; `radixSortIndices` 22.5 and 188 ms). z, not the distance to the eye: moving the camera keeps the order, only
  turning changes it.
- **Where**: `Renderer.update` (per view, before any pass), for the visible splat objects, `updateSplatOrder`
  (`gl/renderable/gaussian-splats-order.ts`): an instance sorted again where its z row's direction turned by more than
  about 1.8° (cosine 0.9995) or its centers or transforms changed, the renderable updated where the order changed. One
  order for all views, no per-view cache: the eyes of stereo look the same way (within the threshold), picks and other
  variants need no order (the nearest fragment above the cutoff wins).
- **Read**: `tSplatOrder`, a data layer (`rgba8`: a `u32` index in four bytes; 4 bytes a splat), per instance a block
  (`aInstance * uVertexCount + slot`, instance ids by the transforms' rows: an instance grid's order kept); the vertex
  shader maps its draw slot `VertexID / 6` through it to the splat, in every variant. On both backends one path, no
  WebGPU storage buffer.
- **Instances**: each instance sorted in its block, drawn one after another (instances overlapping each other composite
  in instance order: accepted); splat objects are never merged (`MergedTypes` unchanged).
- **Cost**: on the CPU, the sort while turning (about every frame: 8-10 ms at 1 M splats, 42-47 at 5 M) and the
  order's upload (4 bytes a splat). On WebGPU the GPU sorts (`gl/compute/splat-sort.ts`, below): 0.4 ms at 1 M, 1.0
  at 5 M to a readback (Dawn, RTX 4070 Ti), 0.1 ms of CPU to record. A worker for the CPU sort not needed for now.
- **GPU sort** (WebGPU, once its kernels compiled; the CPU's until then): per instance a counting sort of the same
  16-bit keys, the depth range from the bounds (no reduction): `splat-keys` (each key, counted into 65536 buckets by
  atomics), the prefix sum of the counts (`createPrefixSum`), `splat-scatter` (each splat to its key's next slot by an
  atomic) into `sSplatOrder`, a storage buffer of the renderable read in place of `tSplatOrder` (`dSplatGpuOrder`).
  Equal keys in any order (the atomics'): a sort per turn, so no flicker while still. The centers copied to a storage
  buffer for the kernels (16 bytes a splat). `Renderer.prepareGaussianSplats` sorts, outside of passes (as
  `prepareCull`, by `DrawPass` and the tracing pass), not `update` (called inside passes: picking, marking). A backend
  switch or a restore releases the old context's order (`releaseForeignSplatOrder`, in `transferTextures`): a foreign
  storage buffer would drop the render object, which the shape representation doesn't build again. The centers and
  order buffers are destroyed with the renderable (`disposeSplatOrder`). A precompile list with splats prepares the
  kernels and the programs of the GPU's order (`dSplatGpuOrder` set where the sort runs, not as recorded): the first
  frame sorts on the GPU and compiles nothing.

### Shading

- **Projection (EWA)**: the 3D covariance projected by the Jacobian of the projection at the center, `Σ' = J W Σ Wᵀ
  Jᵀ`, dilated by 0.3 px² for anti-aliasing; the quad spans the 2D Gaussian's principal axes to `kernelSize` standard
  deviations (default 3, the prototype's 2 cut the falloff at 13.5% and showed circles), clamped to 1024 px. An
  orthographic camera has a constant `J` (the prototype handled perspective only).
- **Fragment**: `exp(-½ dᵀ Σ'⁻¹ d)` times opacity, discarded below `alphaCutoff` (1/255) in every variant, blended as
  the blended transparency (straight alpha over an opaque background, premultiplied over a transparent one); fog, marker
  color, clipping and transparency layers as for points.
- **Params**: `sizeFactor` (scales Σ by its square), `kernelSize`, `alphaCutoff`; the prototype's `splatDebug` (fixed
  circles) as a test helper only, not a param.
- **View-dependent color** (spherical harmonics, degrees 1 to 3): per splat a texel of `rgba8` per basis function past
  the zeroth (`tSplatSh`: rgb the channels' coefficients over the object's largest magnitude, `uSplatShScale`; 12, 32
  or 60 bytes a splat), evaluated in the vertex shader (3D Gaussian splatting's basis and constants) for the direction
  from the camera to the splat in its own frame (the model-view's rotation undone by its transpose), added to the base
  color in the fragment (the color and tracing variants), not below 0. Bytes, not floats (180 bytes a splat at degree
  3): the compressed layout stores bytes too. Turned upright with the splats by the shape's transform (evaluated in
  the splats' frame). The
  degree read is the PLY provider's `shDegree` (0 by default: memory), at most the file's.
- No lighting (splats carry baked radiance: `ignoreDefineUnlit`, as the prototype).

### Passes

- **Draw fact** `GaussianSplats` (bit 28), set by `dGeometryType`, in `OpaqueBlockers` and out of the transparent
  passes (WBOIT, DPOIT, blended transparent, transparent depth); counted as half transparent in `writeContributions`, so
  `opacityAverage` opens the transparent branch.
- **One renderer pass** `gaussianSplats` (`Renderer.renderGaussianSplats`; variant color, the blended transparency's
  states: straight alpha over an opaque background, premultiplied over a transparent one; depth test against the opaque
  depth, no depth write), gated by `Scene.hasGaussianSplats` (the fact in a visible row): after the OIT composite
  (after `wboit.render`, `dpoit.render`, begun on the target with `'opaque-texture'` depth), in blended mode in the
  blended transparency's pass right after it (the drawing buffer's too, on its own depth), in illumination after its
  transparent block on `transparentWithDepth`. Volumes before splats (splats in front of a volume are the common case in
  a mixed scene; volumes inside splat clouds are not). Collected in the prepare plan as any renderer pass.
- **Exclusions**: not occluders (`frame-render.ts`'s occluder pick tests `opaque`, false for splats), out of the
  opaque, tracing, emissive and transparent depth passes, no shadows or SSAO from them (no depth written); outlines and
  DOF ignore them (accepted for a first version). Picked and marked (the pick and marking variants) only where a
  fragment's coverage reaches `pickingAlphaThreshold` (0.5 by default): not by the faint tails of splats in front of
  what is seen; a splat of lower opacity is not picked at all.
- **Transparency mode**: the splat render item is created blended and ignores `setTransparency` (one program, not
  three); the precompile plans them blended too (`Canvas3DPreparer`), its stand-in scenes with `hasGaussianSplats`.

### Formats

- **Layers**: `io` parses, with no or minimal processing and no copies (`reader/ply`, `reader/splat`: the bytes and a
  `DataView`; `reader/spz`: the header and its sections as views, inflated where gzipped); `model` interprets
  (`formats/shape/gaussian-splat/`: `GaussianSplatData`, each splat decoded on demand by `getSplat` from the file's own
  data, the view-dependent color by an accessor; the file's frame, `upright`; the SH constants); `graphics` builds
  (`formats/shape/gaussian-splat.ts`: `buildGaussianSplats` writes each splat once into the layouts the shaders read,
  the shape providers). No intermediate arrays between the file and the geometry (1.5 M splats with degree 3 would be
  over 300 MB of floats).
- **PLY**: `isGaussianSplatPly` (vertex `f_dc_0`, `scale_0`, `rot_0`; or a `chunk` element with `packed_*`) branches in
  `ShapeFromPly` to the splat provider (`model/src/formats/shape/gaussian-splat/ply.ts`): standard (DC color
  `0.5 + C0 · f_dc`, opacity by sigmoid, the log scales as stored, quaternion w-first) and compressed (256 splats a
  chunk, 11-10-11 positions and log scales within the chunk's bounds, smallest-three 2-10-10-10 rotations, 8-8-8-8
  color and opacity, colors within optional bounds), as PlayCanvas's `gsplat-compressed-data.js` decodes it: the
  rotation's index the largest component, 0 for `w`, the first field `w` otherwise (the prototype's mapping was wrong).
  One extension, one provider: no `isApplicable` guessing. `upright` (on by default) turns the scene 180° about x by
  the shape's transform (`getSplatTransforms`, before the file's own transforms; the data kept as stored, each reader
  saying its frame: `GaussianSplatData.upright`, true for SPZ only): 3DGS scenes are y-down, the camera convention of
  their training, upside down as stored (other viewers turn them too). Toggling it builds a new shape of the same
  geometry.
- **Binary tables fast** (`readUniformBinaryTable`): tables of one type in the host's byte order (3DGS: all `float`;
  compressed: `uint` vertices, `float` chunks) read as one typed view (a copy where misaligned), de-interleaved by
  plain loops; the per-value `DataView` loop stays for mixed tables.
- **`.splat` and `.spz`** (`io/src/reader/splat`, `io/src/reader/spz`, `model/src/formats/shape/gaussian-splat/`;
  providers `splat`, `spz` by extension, straight from the bytes to a shape provider, `ShapeFromGaussianSplatFile`; the
  PLY provider's params):
  - `.splat` (antimatter15, no header, 32 bytes a splat): position and linear scales `float32`, color and opacity
    (sigmoid's) bytes, rotation bytes `q * 128 + 128` as `w, x, y, z`; the 3DGS PLY's frame.
  - `.spz` (Niantic) versions 1 to 3 (gzip, read compressed or not: downloads may be inflated already): positions 24-bit
    fixed point (version 1 `float16`), opacity sigmoid bytes, colors bytes of the zeroth harmonic by 0.15, log scales
    `b / 16 - 10`, rotations 3 bytes (version 2) or the smallest three of 10 bits (version 3), harmonics bytes
    `(b - 128) / 128` coefficient-major (degree 4 read to 3); stored right-up-back, the PLY's frame turned 180° about x
    (upright as stored). Version 4 (zstd streams) is refused.
  - `.ksplat` not planned.
- **SOG** (PlayCanvas, step 9, done): a zip (`.sog`) or a folder of `meta.json` and lossless WebP images, a pixel
  a splat in row-major order (`means_l`/`means_u`, `scales`, `quats`, `sh0`, optionally `shN_centroids`/`shN_labels`);
  versions 2 and 3 (`meta.version`), version 1 (no `version`, per-channel ranges) later if needed. The 3DGS PLY's
  frame. Decoding:
  - means: 16 bits an axis (`u << 8 | l`) within `mins`/`maxs`, then `sign(n) * (exp(|n|) - 1)`;
  - scales: log scales from a 256-entry `codebook` by byte;
  - quats: version 3 a 24-bit exponential map (`p = (rgb - 128) / 127`, `theta = π/2 |p|`, `q = (cos theta,
    sin theta p / |p|)`, `w, x, y, z`); version 2 the smallest three (rgb `(b / 255 * 2 - 1) / √2`, alpha 252 to 255
    the largest's index);
  - sh0: the zeroth harmonic from the `codebook` (base color `0.5 + SH_C0 dc`), opacity version 3 `(a / 204)²`,
    version 2 `a / 255`; version 3's opacities above 1 (to 1.5625) clamped to 1 for now (decided; see Open);
  - shN: a palette (k-means centroids, 64 a row, `coeffs` pixels an entry, rgb bytes into the shared `codebook`) and a
    label a splat (24 bits, version 2 16; version 2's `meta.json` gives no `count` nor `bands`: from the palette image's
    size), 1 to 4 bands (read to 3), band order within a channel: the view-dependent
    color's accessor reads the palette (no expansion per splat: the palette is the compression).
  - WebP (decided: the browser's decoder for now): `createImageBitmap` with `premultiplyAlpha: 'none'` and
    `colorSpaceConversion: 'none'`, the pixels read through a WebGL2 texture (`readPixels`), not a 2D canvas (which
    premultiplies: `sh0`'s alpha is the opacity, its rgb indices would be corrupted below 255). The decoder given to
    the reader (`io`, a function of the bytes to RGBA), so tests give theirs; Node (headless) reads no SOG for now. A
    VP8L decoder in TypeScript (the same bytes in the browser and Node) later, if headless SOG is needed.
  - Layers as for the others: `io/src/reader/sog` (the zip or the files, `meta.json` parsed, the images decoded),
    `model/.../gaussian-splat/sog.ts` (the accessor), a `sog` provider (`.sog`; a folder by its `meta.json` later).

### SOG on the GPU (step 10, in progress)

SOG's images are already the layout a shader can read (a pixel a splat, 4 bytes an image); today they make a round trip:
the browser decodes each WebP, `readPixels` brings it to the CPU (`decodeWebpInBrowser`), the model decodes each splat,
graphics encodes it into `tSplat` and `tSplatSh` and uploads them. The view-dependent color is expanded per splat (15
texels, 60 bytes at degree 3) where SOG stores a palette (≤ 64K entries) and a label a splat: at 1 M splats of degree
3 about 100 bytes a splat on the CPU and as much on the GPU, against about 40 on the GPU read in place.

- **Baseline** (today's path, Chrome, RTX 4070 Ti, medians of 3): skull (247K splats, harmonics of degree 3) loads in
  175 ms on WebGL2 (the zip and the WebP decode with readback 62 ms, the build 50 ms), degree 3 adds 300 ms (the
  build 224: the harmonics expanded per splat); 41 bytes a splat on the CPU and 32 on the GPU at degree 0, 101 and 92
  at degree 3; the decoded images 10 MB in passing. Apartment (661K, no harmonics): 346 ms (decode 82, build 146),
  41 bytes a splat on either side. WebGPU's times alike (its decode on WebGL2 contexts of its own, slower on the
  skull: 167 ms); its GPU bytes not measured (the storage pages are not freed with the object), its layout is
  WebGL's. Expected read in place: about 20 bytes a splat on the CPU (positions, order) and 28 on the GPU plus the
  palette, the build and the expansion gone. Large files with harmonics: AlterEyes' swallowtail (3.9 M splats) loads in
  1.6 s on WebGL2 (decode 345 ms, build 885), degree 3 adds 3.0 s (the build 2.3 s) for 398 MB on the CPU and as much
  on the GPU (101 bytes a splat; 99 MB of decoded images in passing); a strawberry (1.5 M): 0.6 s, degree 3 1.4 s more,
  151 MB a side. Read in place the swallowtail would hold about 80 MB on the CPU and 110 MB plus the palette on the GPU,
  its first frame the decode and the means'.
- **PlayCanvas**: the engine reads SOG in place (`gsplat-sog-resource.js`, `gsplat-sog-data.js`, the `sog` vertex
  chunk): each image an `RGBA8` texture from `createImageBitmap` (no premultiplying, no color conversion, no mipmaps),
  read by `texelFetch`; the three codebooks in one `256 × 1` `RGBA32F` texture (scales, sh0, shN); the centers for its
  CPU sort (a worker) by a GPU pass decoding the means into an `RGBA32U` target, read back once
  (`generateCenters`); the bounds from `meta.json` (`mins`/`maxs` mapped out of the log domain); the `ImageBitmap`s
  kept for a device loss. SuperSplat, an editor of splats, decodes SOG on the CPU (`splat-transform`) into its own
  uncompressed textures instead.
- **Here**: the images as textures read in place, a second encoding of the geometry (`dSplatEncoding`: `packed`, today's
  `tSplat`; `sog`). The shader decodes the positions from the means too (`tSogMeansLow`, `tSogMeansHigh`, the
  ranges as uniforms; the group the splat's index): no `tPositionGroup` on the GPU (16 bytes a splat). The CPU still
  needs them for the CPU sort, loci bounds, the position iterator and a switch to WebGL: the geometry's
  `positionGroup` from the means alone (two images read back and decoded as today, 8 bytes a splat: no float render
  target, which WebGL1 and some WebGL2 lack, and the decode tested), when depends on the backend:
  - WebGL: at load, before the first frame (its sort is the CPU's);
  - WebGPU: in the background after the first frame, off the load's path. The GPU sort's centers buffer is filled by
    a compute pass from the means textures (`splat-centers`, with the sort's kernels in the precompile), not from
    `positionGroup`. Until the means arrive, loci bounds are the object's (`meta.json`) and the position iterator
    has none; until the sort's kernels compile (no precompile list), the splats draw in the order stored, not on the
    CPU.

  Loci bounds pad by the largest scale of the codebook (no per-splat scales on the CPU; looser than `tSplat`'s).
  Everything else is not read back:
  - values: `tSogScales`, `tSogQuats`, `tSogSh0`, `tSogShLabels`, `tSogShCentroids` (`image-uint8` textures from
    `ImageBitmap`s, their own width: `i % w, i / w`), `tSogCodebook` (`256 × 1` `rgba32f`: scales, sh0, shN, as
    PlayCanvas), `dSogVersion` (2 or 3: rotation and opacity), the palette's coefficients an entry and count, the
    degree (`uSplatShDegree`, at most the file's);
  - shader: one decode function per encoding giving rotation, scales, opacity, base color and the coefficients
    (from the palette by label: no expansion); the projection, order, picking and marking unchanged. GLSL (WebGL1:
    `texture2D` by float index math, exact to 2^24) and WGSL (`textureLoad`);
  - color: the base color is the file's, read in the shader: a color theme of the file's colors with the `direct`
    color type (as direct volumes), the other themes and overpaint as today (`tColor`); the opacity from `sh0`'s alpha;
  - upright (decided): the shape's transform, for every format (below), so the sog decode has no turn;
  - bounds from `meta.json` (exact, before the means are read), the group bounds from `positionGroup` once read.
- **Textures from `ImageBitmap`s**: WebGPU's textures load them (`copyExternalImageToTexture`), WebGL's need it
  (`texImage2D` with the unpack state of `decodeWebpInBrowser`). The bitmaps closed once uploaded. A context restore
  or a backend switch decodes again from the file's bytes, which the state tree holds anyway (decided: nothing kept,
  unlike PlayCanvas's bitmaps at 4 bytes a pixel each; the restore asynchronous, the representation built again for
  the new context as `GeometryUtils.isCompatible` does, the splats missing for its frames).
- **Upright as a transform** (done, every format; see PLY above): the per-format turn is gone from the readers, the
  harmonics stay right (evaluated in the object's frame, the view direction turned back by the model-view's
  transpose; checked against the images of the baked turn: the cactus at degree 3 from two sides alike).
- **Premultiplying**: the SOG spec warns that some browsers (WebKit under WebGPU) premultiply despite
  `premultiplyAlpha: 'none'`, corrupting `sh0`'s rgb indices below alpha 255. A probe at load (a 1 × 1 WebP of alpha
  128 through the same path, its rgb compared) falls back to today's CPU path, which keeps the CPU reader alive (also
  for headless reading with a decoder given, and the tests).
- **WebGPU bindings**: five to six more sampled textures in the splat programs; their count against
  `maxSampledTexturesPerShaderStage` (16 by default) with the global and color textures checked first. The codebook as
  a uniform array is the alternative (4 KB: a block on WebGL2 and WebGPU, too large for WebGL1's classic uniforms).


## Salvage from the v5 prototype

| Prototype                                            | Here                                                                           |
| ---------------------------------------------------- | ------------------------------------------------------------------------------ |
| `gaussian-splat-ply.ts` (standard and compressed)    | ported, writing the arrays directly; the compressed rotation's mapping fixed   |
| `computeCovariance` and its tests                    | `getSplatCovariance`, the shaders' oracle (bounds from the largest scale)    |
| EWA math in `gaussian-splats.vert.ts`                | ported to GLSL (v6 chunks) and WGSL; orthographic added, premultiplied output  |
| `gaussian-splats-sort.ts` (16-bit counting sort)     | kept (faster than `radixSortIndices`), keyed by view z; its tests rewritten    |
| own textures (`tCenterGroup`, `tCovariance0/1`, ...) | dropped: `tPositionGroup`, `tSplat`, `tSplatOrder` data layers                 |
| sort in `Renderer.renderObject`, `WeakMap` cache     | moved to `Renderer.update` (once per view, before the passes)                  |
| `renderBlendedGaussianSplats` after the OIT resolve  | kept as the idea: a `PassDraws` row and a sub-pass in each mode                |
| v5 binary PLY parser (other stash)                   | not needed: v6's parser reads binary                                           |

## Steps

1. **Geometry and data** (done): `geo/geometry/gaussian-splats/` (the geometry's arrays, bounds padded by 3σ of the
   largest scale, group mapping), the `tSplat` encoding with its round-trip and covariance tests, the page bits by the
   schema's layers (`webgpu-data-layers.test.ts`). The kind's wiring moved to step 2: `GeometryUtils`,
   `RenderObjectValues` and `RenderablePrograms` need the kind's shader code once it is in the unions, and a shader code
   its programs (the enumerating tests).
2. **Kind and shaders, unsorted** (done): `GaussianSplats.Utils` (values, position iterator from `tPositionGroup`,
   params), the schema (`gl/renderable/gaussian-splats.ts`, `dGeometryType`), the kind in `geometry.ts` and
   `render-object.ts`; `gaussian-splats` programs in GLSL and WGSL (`QuadItems`, all variants), registered on both
   backends and in the enumerating tests (`graphics-shaders.ts`, `shader-parity`); the draw fact and the
   `gaussianSplats` pass in each transparency mode and illumination (`draw-table.test.ts`: the routing). A headless
   WebGPU test (`webgpu-gaussian-splats.test.ts`): a splat's profiles against the projected Gaussian (isotropic;
   stretched and turned).
3. **Order on the CPU** (done): view z keys, a counting sort, `tSplatOrder` on both backends, sorted in
   `Renderer.update` where the view turned (`gaussian-splats-order.test.ts`: the order for a camera, kept moving, again
   turned, per instance; `webgpu-gaussian-splats.test.ts`: two splats from either side); measured 1 M and 5 M splats,
   the worker deferred.
4. **PLY** (done): detection, both layouts, the typed-array path in the parser; the provider in `ShapeFromPly`; tests
   with inline fixtures (a few splats of each layout, as `ply.test.ts` builds its data; the compressed one by an encoder
   after PlayCanvas's decoder). Two public compressed scenes loaded in the Viewer (WebGL2, WebGPU); the standard layout
   only by the tests, Prism not tried.
5. **Picking, marking, themes** (done): pick and marking above the picking threshold, loci bounds by group
   (`geo/shape/shape.ts`, padded by 3σ of the group's largest scale); tests: a pick by coverage on WebGPU, the group
   bounds, visual cases `highlight-select` and `overpaint-transparency` (WebGL2 and WebGPU alike), `sphere` in the
   smoke subset.
6. **GPU sort (WebGPU)** (done): a counting sort by atomics on the prefix sum, into the storage buffer the vertex
   shader reads, no readback; the CPU order kept for WebGL and until the kernels compiled. Tests: 200k splats in two
   instances against their depths (`webgpu-splat-sort.test.ts`), the splats in order from either side drawn with the
   GPU's order; measured 1 M and 5 M splats.
7. **Spherical harmonics** (done): degrees 1 to 3 as a layer of bytes, evaluated in the vertex shader for the view
   direction; `shDegree` (0 by default). Tests: both layouts' coefficients (channel-major; compressed `sh` bytes) and
   their turn, a splat's red by the direction on WebGPU, a visual case (WebGL2 and WebGPU alike); by eye on the cactus
   and strawberry (highlights).
8. **More formats** (done): `.splat`, `.spz` (versions 1 to 3). Tests: files written by the reference packers'
   layouts (`convert.py`, `load-spz.cc`) read back, gzipped or not, upright or not, `.spz` harmonics; by eye on
   antimatter15's nike (`.splat`) and Niantic's horned lizard (`.spz`).

9. **SOG** (done): `io/src/reader/sog` (zip, `meta.json`, versions 2 and 3, the images by the given decoder,
   `decodeWebpInBrowser` by default), `model/.../gaussian-splat/sog.ts`, the `sog` provider. Tests: zips of fake images
   (the decoder given), splats encoded into images as splat-transform's writer does (both rotations, the clamped
   opacity, the palette by label); by eye on PlayCanvas's skull (version 2, with harmonics: the bronze's highlights) and
   apartment, WebGL2 and WebGPU.

10. **SOG on the GPU** (in progress, see the section; whole, not the palette in the packed encoding first: one new
    encoding, no intermediate layout). Done:
    - the baseline (a) and upright as the shape's transform for every format (d);
    - textures load `ImageBitmap`s on both backends, not kept as their source (closed once uploaded);
    - the `sog` encoding (`GaussianSplats.createSog`, `buildSogSplats`): SOG's archive read without decoding
      (`parseSogArchive`), the means alone decoded and read back (`decodeSogImage`: the CPU's centers), the other
      images uploaded as the browser decodes them; the decode in GLSL (WebGL1 by coordinates) and WGSL, positions from
      the means (`tPositionGroup` a placeholder, the sorts' centers `splatCenters`), the palette read by label (e);
    - the file's colors as the `direct` color type (`ShapeUtils.getTheme`; the splat shaders take `vColor` for it);
    - the textures the geometry's: destroyed when replaced or with the representation (`ShapeProvider.dispose`), a
      new context building them again through the getter (`ShapeGetter` given the context; the shape representation
      calls it again where `isCompatible` fails and there is no `makeCompatible`); a representation that is not the
      geometry's (a provider's other one, e.g. recreated by `TransformShapeProvider`, which passes the context on)
      builds its own (f).

    Measured (Chrome, WebGL2, medians of 3): the swallowtail loads in 670 ms (1563 before), degree 3 adds 561 ms
    (3042); the strawberry 323 ms (596); 21 bytes a splat on the CPU (101 at degree 3 before), on the GPU 30 to 40 at
    degree 3 (the palette included; about 100 before). The images as the packed path's within the packed layout's
    quantization (skull: 0.04% of pixels past 8 levels), WebGL1 as WebGL2, a backend switch drawn again.
    Tests: a synthetic SOG drawn in place as its CPU decode in the packed layout (headless WebGPU, within 8 levels);
    the geometry's textures owned and released; the archive read without decoding; the getter given the context;
    the transform passing it on.

    Then:
    - WebGPU reads no image back at load: the GPU sort's centers decoded from the means textures (`splat-centers`, one
      of `SplatSortKernels`, so precompiled with them), the CPU's read after the first frame (`setSogCenters`); until
      then the bounds are the means' ranges' (`centersBounds`) and a splat's loci bounds the whole's. The swallowtail's
      first frame 558 ms on WebGPU (759 reading them at load). The camera's first fit is to those bounds (the box of
      the ranges: up to √3 the centers' sphere; the skull's 3.4 against 2.8), so it frames a little wider than WebGL
      (decided: accepted, not re-fitted when the centers arrive).
    - the premultiplying probe (a 1 × 1 WebP of alpha 128 uploaded and read back, once per context, while the file is
      read), falling back to the CPU path with a warning;
    - the image size the means texture's own (a bitmap's is 0 once closed: the shaders read 0 × 0 images otherwise);
    - tests: the precompile of the `sog` programs (headless WebGPU), the GPU's centers as the CPU's, the centers read
      later (groups mapped, the bounds the ranges', a splat's loci bounds the whole's until they arrive); the visual
      case `geometry/gaussian-splats/sog` (WebGL1, WebGL2, WebGPU: a sphere encoded as version 3 with a palette).

    Left: the browser's probe failing is not exercised (Chrome does not premultiply); the deferred path has no visual
    case (it needs SOG files as stored: WebP encoded, the browser's decoder).

## Open

- SOG version 3 opacities above 1 (to 1.5625: `(a / 204)²`) are clamped to 1: `tSplat`'s opacity is unorm16, the
  `sog` encoding's shaders clamp as the CPU's decode. Keeping them needs a wider opacity (e.g. a scale per object) and
  the blending of an opacity above 1 (as PlayCanvas does).

- Splats and other transparent objects in one scene composite as "all other transparency, then splats": exact only
  where they don't interleave in depth. Sorting splats into WBOIT/DPOIT is possible (DPOIT peels them correctly at its
  cost per layer) but not planned.
- Depth for post-processing (proposed, not decided): splats write no depth, so SSAO, outlines and DOF ignore them.
  (A, preferred) splats in the transparent depth pass (`depthTransparent`, the depth SSAO and outlines read for
  transparent objects): a splat's center depth where its coverage reaches the picking threshold, as picking; no sort,
  no target, drawn only where post-processing needs that depth; a flat disc a splat (steps on large splats in SSAO).
  (B) an expected depth: a float target in the splat pass blended as the color, resolved to depth where the opacity
  passes a threshold: smooth, but every splat pipeline and a pass more. Whether DOF reads the transparent depth is to
  check.
