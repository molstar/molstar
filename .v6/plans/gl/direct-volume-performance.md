# Plan: direct-volume rendering performance

## Status

- **Current**: steps 0-5, 8, 10 and 11 done (the Gaussian volume keeps `jumpLength`, see step 4); steps 6, 7 (the data
  layout) and 9 open, not worth it on the numbers below.
- **Measurements** (RTX 4070 Ti, Chrome, 2560x1440, `node scripts/benchmark.mjs --filter 'render/direct-volume/*'
  --gpu-timing --runs 7`; the cases in `tests/src/benchmark/direct-volume.ts`, on demand). GPU time of the pass the
  volume draws in (`wboit transparent`, the default transparency), medians of 7 runs of 30 frames, WebGL2 / WebGPU:

  | case                         | baseline     | step 1       | step 2       | step 4       | step 3       | step 5       | final |
  | ---------------------------- | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | --- |
  | `blobs-256-dense`            | 6.4 / 8.3 ms | 3.7 / 4.3 ms | 3.5 / 3.3 ms | 3.3 / 3.0 ms | 2.4 / 2.3 ms | 0.9 / 0.8 ms | 0.8 / 0.9 ms |
  | `map-256`                    | 8.5 / 8.5 ms | 3.8 / 4.4 ms | 3.0 / 2.9 ms | 1.0 / 1.1 ms | 0.8 / 0.8 ms | 0.4 / 0.4 ms | 0.4 / 0.5 ms |
  | `map-256-zoom`               | 23.3 / 33.8  | 16.2 / 19.5  | 11.8 / 12.8  | 4.0 / 4.2 ms | 3.2 / 3.4 ms | 1.1 / 1.1 ms | 1.1 / 1.2 ms |
  | `map-256-float`              | 15.9 / 15.0  | 13.9 / 14.4  | 13.0 / 13.5  | 1.6 / 1.6 ms | 1.1 / 1.1 ms | 0.5 / 0.5 ms | 0.5 / 0.6 ms |
  | `map-256-highlight`          |              | 4.8 / 5.2 ms | 3.4 / 3.0 ms | 1.1 / 1.2 ms | 0.9 / 0.9 ms | 0.4 / 0.5 ms | 0.4 / 0.5 ms |
  | `tomogram-512x512x96`        | 8.9 / 8.7 ms | 2.2 / 2.4 ms | 1.8 / 1.8 ms | 1.2 / 1.3 ms | 0.9 / 1.0 ms | 0.4 / 0.5 ms | 0.4 / 0.5 ms |
  | `map-256x256x64-anisotropic` | 7.6 / 7.9 ms | 3.9 / 5.5 ms | 3.0 / 3.0 ms | 2.1 / 2.2 ms | 1.5 / 1.6 ms | 0.6 / 0.7 ms | 0.6 / 0.7 ms |
  | `7qpd-gaussian-volume`       | 8.0 / 9.0 ms | 4.6 / 4.8 ms | 3.2 / 3.4 ms | 3.8 / 3.3 ms | 2.8 / 2.4 ms | 1.3 / 1.2 ms | 1.0 / 1.2 ms |

  - Baseline: everything else in the frame is below 0.5 ms. `float` costs about twice `byte` (finding 8 matters more
    than first thought, and after step 1 it dominates `map-256-float`). The Gaussian volume runs with `jumpLength` 4
    (its default).
  - A first baseline was 2-4x higher: the benchmark's theme lacked `locationKinds`, so `volume-value` colored per cell
    (`vertex`, eight texel reads per sample through `texture3dFrom1dTrilinear`) instead of by value in the shader as in
    the plugin (`direct-location`). The test helper `volumeRepr` now passes the provider's `locationKinds` (finding 10).
  - Step 1: images identical to the baseline (at most 2 levels), but for 21 pixels of the anisotropic case where rays
    now reach the volume's end. The steps before the box were most of the cost of the tilted slab (-75 %); the
    zoomed-in case gains least (-30 % / -42 %): its rays start inside or near the box.
  - Step 2: images identical to step 1 (0 levels), but for the highlighted case, where the haze over the whole box is
    gone (`map-256-highlight` added with step 2: markers colored in the volume shader need the marking pass off).
    -10 to -30 %, more on WebGPU, which was 5-40 % slower than WebGL2 after step 1 and is on par now.
  - Step 4 (volumes): images identical to step 2 (0 levels). Sparse maps 3x faster, the float map 8x (empty space is
    no longer read: bandwidth was its cost); the dense blobs and the Gaussian volume (no bricks yet) unchanged. The
    brick ranges cost about a quarter of the grid texture's build (256³: 204 ms against 764 ms, 512x512x96: 320 ms
    against 1399 ms). The Gaussian volume without `jumpLength` (its default 4) costs 5.6 / 7.4 ms instead of
    3.8 / 3.3 ms; with bricks of its own it is still slower than with `jumpLength` (step 4).
  - The Gaussian volume's times vary by up to 35 % between page loads (3.7 and 5.0 ms for the same build on WebGL2), the
    volumes' by a few percent: compare it within one page, and over several loads.
  - Step 3: jitter alone changes images by at most 3 levels at 3 steps per cell; the default of 2 (was 3) is
    indistinguishable from 3 and 10 in all cases (also the Gaussian volume zoomed 3x), -25 %; at 1 the jitter shows as a
    fine grain on the Gaussian volume (32 k pixels off by more than 8 against 10), multi-sampling halves it.
  - Step 5 (the column: the volume at half resolution, its WBOIT accumulation, evaluate and composite): the whole frame
    34-63 % faster than at full resolution in the same page (`map-256-zoom` 3.5 -> 1.4 ms, both backends), the rays
    3.1x cheaper (the ideal at a quarter of the pixels 4x), the composite 0.04 ms. Against full resolution the maps differ
    by at most 9 levels (0-3 pixels off by more than 8), the Gaussian volume by up to 24 (fine density), the volume at
    the edges of an opaque isosurface by up to 71 in single pixels along the edges (`map-256-isosurface`, added with
    step 5); side by side the edges look the same. DPOIT and blended mode and stereo checked by hand (images).
  - `blobs-256-dense` shows ring banding: not from the steps (jitter doesn't change it) but from the 8-bit values. Its transfer alpha (0.01, stored as 2/255 in the
    8-bit transfer texture) is below the lighting threshold: the case measures unlit samples.
  - `multiSample` variants not added: the volume's cost multiplies by the sample count (the scene is drawn per
    sample), nothing to learn from measuring it.
  - Single-run noise: `map-256-zoom` on WebGL2 varies by up to 40 % between runs; compare medians over several runs.
  - Final (after steps 8, 10, 11; the step 5 column's passes): level with step 5 but for the lit cases on WebGL2 (the
    normal matrix per instance, step 8: the Gaussian volume 1.3 -> 1.0 ms). The whole frame with a volume 0.7-1.5 ms
    on both backends (from 6.5-34 ms at the baseline); with the marking pass on and the volume highlighted 1.3 / 1.4 ms
    (`map-256-highlight-marking`, step 11). Full visual suite (all profiles), the quick, medium and heavy lanes pass
    (two heavy tests of `packages/sequences` fail on this machine, unrelated: a symlink needs rights Windows withholds,
    and the unused-exports listing).
- Decided: no marker glow over (almost) fully transparent voxels; `jumpLength` goes from volumes, which skip bricks,
  and stays with the Gaussian volume, where bricks measured slower (step 4); a slightly different look (lower default
  `stepsPerCell`, jitter) is acceptable.
- Scope: the ray marcher (`gl/webgl/shader/direct-volume.frag.ts`, `gl/webgpu/shader/direct-volume.wgsl.ts`, the two
  are line for line the same loop), its geometry (`geo/geometry/direct-volume/`), the volume representation
  (`repr/volume/direct-volume.ts`) and the Gaussian volume (`repr/structure/visual/gaussian-density-volume.ts`, packed
  groups), and the passes that draw it (`canvas3d/passes/draw.ts`, `illumination.ts`, `multi-sample.ts`).

## Where things stand

- **Draw**: one box per instance, culling `none`, front faces discarded in the shader, `gl_Position.z` at the far plane.
  Below `volumeResolutionScale` 1 (Canvas3D, default 0.5) the volumes go through `VolumePass`
  (`canvas3d/passes/volume.ts`, step 5) instead of the passes below. The `volume` pass (premultiplied blending), or a
  WBOIT layer of their own in `wboit` mode, under the transparent objects where there are any (step 10), else over the
  postprocessed image. The `color`, `marking`, `pick` and `depth` variants (step 11; not `emissive`, `tracing`).
- **Ray setup per fragment**: `inverse4(vTransform)`, the normal matrix per instance from the vertex shader (step 8);
  the samples at
  `start + t * step` (start at `max(near, bounding-sphere near)`) in unit grid coordinates and view space, affine in
  `t`; `t` from the box entry (slab test) to the box exit, the far plane or the opaque depth (as view z), whichever is
  first; the loop also ends at alpha 0.95 (step 1).
- **Per step**: the brick test (step 4), the grid fetch (trilinear, RGBA: normal in RGB, value in A), the transfer
  function fetch; then, only where the transfer alpha is above 0 (step 2), the group index, the color (texture types: a
  fetch per step, `vertex` types eight), overpaint, marker fetch, fog, blending;
- **Lighting** only where the transfer alpha is at least 0.01; packed groups (Gaussian volume) compute the gradient
  from six more grid fetches.
- **Empty space**: volumes skip bricks the transfer function makes transparent (step 4); the Gaussian volume
  probes `jumpLength` (4 Å) ahead where the density is below 0.01 (approximate: may step over thin features).
- **Steps**: `stepsPerCell` 2 by default (step 3), step `min(cellDim) / stepsPerCell`; `uMaxSteps` the steps along the
  box's longest diagonal (a loop bound only, the box exit ends the ray); the samples offset per pixel (jitter, step
  3); opacity per step corrected exactly (finding 14).
- **Data**: RGBA per voxel: 4 B (`byte`), 8 B (`halfFloat`), 16 B (`float`, 2 GiB at 512³); normals precomputed on the
  CPU (`createVolumeTexture3d`, read at `Tensor.strides`, step 7) on every geometry build, with the brick ranges;
  transfer function changes upload 256 bytes and the brick visibility.
- **Multi-sampling**: `multiSample` re-renders the whole scene per sample, the volume included: its cost times the
  sample count while the view is still.
- **Coverage**: visual cases `geometry/direct-volume/` (`blobs` at half and `blobs-full-resolution` at full resolution,
  `opaque-edges`: an opaque isosurface inside, `bricks`: a narrow transfer function, `gaussian-volume`) and
  `transparency/{blended,wboit,dpoit,dpoit-4}/volume` (an opaque box inside, per mode); the
  benchmark cases `render/direct-volume/*`. The ordering against transparent objects (finding 13):
  `transparency/*/volume-behind*` (step 10); marking and picking: `renderer/marking/direct-volume`,
  `renderer/pick/direct-volume` (step 11).

## Findings

1. **Steps wasted before the box** (largest suspected cost for elongated or rotated grids): the ray starts at the
   bounding sphere and each step until the box costs a `mat4` transform. A tomogram of 1000x1000x200 cells has a radius
   of about 720 cells against a half-depth of 100, up to ~1800 empty steps per ray at 3 steps per cell.
2. **Rays truncated** (correctness, same cause): `uMaxSteps` covers the grid diagonal in cells, but the step is the
   smallest cell edge, and the steps before the box count too. Anisotropic cells (e.g. 1x1x4 Å) or a far sphere-near
   start end rays early, cutting off the far part of the volume.
3. **Work on transparent samples**: with the default ramp (from mean + σ) most samples have transfer alpha 0, yet each
   runs the depth test, the group index, color and marker fetches, fog and blending. An early `continue` after the
   transfer lookup skips them. Markers lift alpha to `strength * 0.002`, so highlighted or dimmed volumes glow where
   the transfer alpha is 0; that glow is not wanted and goes with the skip.
4. **Per-step matrix work**: unit position and view position are affine in the ray parameter; marching in unit space
   (start and step transformed once) and keeping the view z as `z0 + i * dz` removes two `mat4` products per step.
   The opaque depth turns into a ray end once (`tMax = min(tExit, tOpaque)`), dropping the per-step depth test and the
   WebGL depth read.
5. **No empty-space skipping**: cryo-EM maps and tomograms are mostly below the ramp. Bricks with their value range
   (min/max per 8³ cells, on the CPU while the texture is built) and a visibility per brick from the transfer function
   (does it reach alpha > 0 in the brick's range; a prefix table over the 256 entries, recomputed on transfer changes,
   ~32 k bricks at 256³) let the ray skip whole bricks.
6. **Steps per cell are high for quality**: banding without jitter and narrow transfer peaks (the sharp-surface
   presets are ~0.012 wide) need 3 or more steps per cell. A ray offset per pixel (interleaved gradient noise,
   varied per multi-sample sample) and pre-integrated transfer functions (a 256x256 table of the integral between two
   values) would allow 1 to 1.5 steps per cell at equal quality, with multi-sampling accumulating the noise away.
7. **Full resolution always**: rays per pixel at the drawing buffer's resolution, re-run per multi-sample sample.
   Volumes are smooth: half resolution with depth-aware upsampling (as SSAO, `ssao-performance.md`) is 4x fewer rays.
8. **Bandwidth**: `float` volumes fetch 16 B per sample, 4x `byte`. A value-only texture (`r16f`/`r32f`, `r8`) with
   normals from central differences only on lit samples (finding 3 makes those rare) or normals in a separate RGBA8
   texture trades fetches for cache footprint; to measure, not to assume.
9. **Minor** (step 8): `inverse4`/`adjoint` per fragment (per instance in the vertex shader or on the CPU instead); the
   `mat4` varying is interpolated on WebGL (`flat` in GLSL ES 3.0); culling front faces by state (with the reflection
   facts, as other objects) instead of `discard`.
10. **Per-cell colors are expensive** (found measuring): a color theme of `vertex` granularity reads eight texels per
    sample (`texture3dFrom1dTrilinear`, manual trilinear filtering from the 2D color texture); 2-4x the cost of the
    palette lookup by value. Themes for direct volumes are by value (`direct`) in the plugin, so only custom themes
    pay it; a 3D color texture with hardware filtering would make it one fetch. Not planned until needed.
11. **NaN pixels** (found measuring, fixed with step 1): `texture3dFrom1dTrilinear` divided `0 / 0` for samples on
    grid coordinates (`ceil == floor`), black pixels in rings and lines; the fraction is now `pos * gridDim - floor`.
    A zero gradient (byte normals filtered halfway between 127 and 128) faces the camera instead of normalizing to
    NaN.
12. **Not highlighted by default, never picked** (found measuring; fixed in step 11): with the marking pass
    on (the default), markers are not colored in the object shaders (`isMaterialColorMarker`), and the marking pass
    draws `scene.primitives` only, never `scene.volumes`, and direct volumes have the `color` variant only
    (`getObjectState`; WGSL `variants: ['color']`): marked direct volumes show no highlight. Without a `pick` variant
    they are never picked either: hovering a volume marks nothing (`getDirectVolumeLoci` exists but is never reached),
    marks come only from loci of other representations (an isosurface or slice of the same volume, whole-volume
    loci). With the marking pass off the volume shader colors its markers itself (step 2).
13. **Transparent objects in front of volumes are drawn behind them** (found reviewing, step 10: volumes now drawn
    under the transparent objects; not from step 5): every
    mode draws the opaque objects, resolves the transparency (blended, WBOIT, DPOIT), postprocesses, and composites
    the volumes last over that (`draw.ts`); the rays stop at the opaque depth only. A transparent surface in front of
    a volume ends up behind it. Illumination does the same in its blended and DPOIT modes; in its WBOIT mode the volumes
    share the transparent objects' accumulation, order independent but weighted by the ray's entry depth (a volume is
    one layer there). Checked with a map, a transparent box in front and one behind (alpha 0.5) in all three modes, at
    both volume resolutions and with illumination: opaque objects, objects behind and several volumes are right.
14. **Opacity per step corrected linearly** (fixed): a sample's opacity was the transfer function's times
    `1 / stepsPerCell` (`uTransferScale`), close to the exact `1 - (1 - a)^(1 / stepsPerCell)` only for small opacities
    (a dense transfer function, 0.6, looked lighter the fewer the steps). Now exact in both shaders (a `pow` per visible
    sample). A transfer function peaking at 0.6, against 10 steps per cell: pixels off by more than 8 at 3 / 2 / 1
    steps 184 / 3148 / 23406 before, 189 / 311 / 4245 after, the mean bias gone; the default cases (opacities up to
    0.05) unchanged (at most 2 levels).
15. **DPOIT textures unbound before the first peel** (found adding visual cases, fixed; not of the volumes): under the
    full visual suite's load `transparency/dpoit-4/volume` failed on WebGL2 ('invalid operation'): a mesh's draw read
    its `tDpoitDepth` unit holding the opaque depth texture, the target's depth attachment. The renderer's shared
    global textures started with `tDepth` and `tDepthAlpha` only, the dpoit ones were added by the first dpoit pass,
    but the `dpoit` mode's programs sample them in every color draw, also in the opaque pass before: their units held
    whatever was bound last (here the volume composite's depth). Now every global texture is in the list from the start
    (placeholders until a dpoit pass sets them); 5.x has the same list. The case is in the suite again (three runs of
    the failing combination pass, the full suite passes).
16. **Holes where rays turn opaque** (found in use, fixed): in the `wboit` mode the transparent pass discards fragments
    of opacity exactly 1, which the opaque pass draws for other objects; direct volumes have no opaque pass, so their
    fully opaque rays left holes (white speckles over a white background). Rare before finding 14's exact correction
    (float rounding of many samples, also in 5.x), a pattern following the jitter after it (a transfer function at 1
    makes one sample opaque). The WBOIT write (both backends) keeps direct volumes' opaque fragments. Visual case
    `geometry/direct-volume/opaque-transfer`.

## Steps

0. **Measure first**. Benchmark render cases (`pnpm benchmark`, `render/*`) at 2560x1440 on WebGL2 and WebGPU, GPU time
   of the `volume` pass: a dense 256³ blob map (most samples visible), a sparse 256³ map (default ramp, ~5 % visible),
   an elongated tomogram-like grid (512x512x96, rotated 30°), an anisotropic grid (cells 1x1x4), a Gaussian volume
   (packed groups). Done: the baseline is in Status; `--images <dir>` writes each case's image, the reference for the
   quality checks of later steps.
1. **Box clipping and a bounded loop** (findings 1, 2, 4): slab test in unit space for `tEnter`/`tExit`, the opaque
   depth as `tOpaque`, march in unit space from `max(tEnter, tNear)` to `min(tExit, tOpaque)`, view z incremental,
   `uMaxSteps` from the box diagonal in step units (cell size aware). Same image, except truncated rays now finish.
   Both backends together; the visual cases keep passing. Done: results in Status; visual cases `geometry/**`,
   `gpu-paths/**`, `fallback/**` pass on WebGL2, WebGL1 (the 2D atlas path) and WebGPU.
2. **Skip transparent samples** (finding 3): `continue` after the transfer lookup when alpha < 1/255, before the
   group, color, marker, fog and blending work; markers no longer lift transparent voxels (the
   `max(alpha, strength * 0.002)` in the volume's marker blending goes). Done, with alpha `<= 0` instead of `< 1/255`:
   the filtered 8-bit transfer function has alphas below 1/255 only at ramp edges, skipping exact zeros keeps the image
   identical. The marker lift went from the shared marker code of both backends (it only served direct volumes). Full
   visual suite passes.
3. **Ray jitter** (finding 6, first half): a per-pixel start offset of `[0, 1)` steps, varied per multi-sample
   sample; then check the quality at lower `stepsPerCell` and lower the default to the lowest count that matches 3
   without jitter. Done: `rayJitter` in both shaders, interleaved gradient noise
   per pixel plus a hash of the projection's off-center terms (the multi-sample sub-pixel offset, 0 otherwise: a
   static pattern without multi-sampling), the samples at `k + jitter` (the brick skip keeps the offset); `stepsPerCell`
   2 by default. Measured against 10 steps per cell (`steps-comparison`, temporary cases, not kept): see Status.
4. **Empty-space skipping** (finding 5): brick min/max texture (`rg8` per 8³ cells, built with the grid texture),
   brick visibility recomputed on transfer changes (CPU, small), the loop advances to the brick's exit when the brick
   is invisible. `jumpLength` removed (`breaking-v6-changes.md`; saved states with it ignore it). The Gaussian volume
   builds its texture on the GPU: its bricks from a reduction kernel on WebGPU (see Compute), without skipping on
   WebGL at first. Done for volumes: `createVolumeBrickRange` (CPU, `repr/volume/util.ts`) per 8³ cells and the next
   cell along each axis, as 0-255 rounded outwards (`BrickRange` on the geometry); the visibility per brick
   (`tBrickTex`, a 2D data texture read with `readFromTexture`, all backends) recomputed in `updateValues` from the
   transfer function; the loop jumps to the first step past an invisible brick's exit face. `jumpLength` is a param
   of the Gaussian volume only (its visuals' params, default 4; `breaking-v6-changes.md`), the geometry reads it where
   given.
   - Gaussian volume bricks built and dropped: ranges rendered from its density on the GPU (a compute kernel on WebGPU,
     whose pipeline is created synchronously, as a render pipeline compiles in the background and skips a one-off draw;
     a fragment pass on WebGL2; none on WebGL1, where the pass wrote nothing without an error), the ranges of volumes
     uploaded into the same texture, the visibility tested in the shader from a 256x256 table of (lo, hi) built from
     the transfer function. Images identical, but slower: the Gaussian density is above 0 near every atom, so nearly
     all of its bricks inside the molecule are visible; bricks alone measured level to 30 % slower than `jumpLength`,
     bricks and `jumpLength` together 25 % slower (4.7 / 4.1 ms against 3.8 / 3.3 ms), and volumes 5-8 % slower with
     the second, dependent read per brick (range, then table).
   - Also measured and dropped: testing a brick once (at its exit) instead of comparing the brick index per step, 3-15 %
     slower; 4³ bricks, volumes up to 2x slower; testing the brick only after a transparent sample, volumes 25-40 %
     slower.
5. **Reduced resolution** (finding 7): the volume pass into a target at `resolutionScale` (0.5 default to measure)
   with depth-aware upsampling into the color target; full resolution when the view is still and multi-sampling is
   on, if the half-resolution image differs visibly. Done: `VolumePass` draws the volumes into a target of `volumeResolutionScale` (Canvas3D and
   image props, default 0.5) as the transparency mode does at full size (weighted blended into WBOIT targets of that
   size in the `wboit` mode, evaluated over transparent black: a premultiplied layer; the volume pass otherwise), then
   composites it (`volume-composite`, both backends): the four nearest texels weighted bilinearly and by how close the
   opaque depth their rays stopped at is to the pixel's (read at the full-size pixel each ray read: the shaders read the
   depth at the fragment's position over `uDrawingBufferSize`, set to the target's size; WGSL scales its integer load).
   Full size, illumination and packed depth (no depth textures) keep the passes as before. The full resolution while
   still was not needed: the half-resolution images differ too little. Prism sets it in the Canvas card's Advanced
   section (`canvas.settings`, Volume resolution).
6. **Pre-integrated transfer function** (finding 6, second half): 256x256 table built on transfer changes; the loop
   looks up `(prevValue, value)` (the unused `prevValue` today), opacity corrected by step length. Only if step 3 does
   not already reach 1-1.5 steps per cell at equal quality.
7. **Data layout** (finding 8): measure value-only plus on-the-fly gradients against RGBA for `byte` and `float` at
   256³ and 512³; change the default storage only with a measured gain. Separately the CPU texture build (perf-audit
   V1): done, `createVolumeTexture3d` and `createVolumeBrickRange` read at `Tensor.strides` (a copy read through the
   accessors for the spaces without: periodic and wrapped volumes) instead of `dataOffset` per cell, the normals inline:
   the texture 2.2-3.5x faster (256³ byte: 0.9 -> 0.4 s, float 0.6 -> 0.23 s; 512x512x96 byte: 1.6 -> 0.65 s), output
   identical (unit test against the build cell by cell, every axis order, a space without strides).
8. **Minor items** (finding 9). Done: the normal matrix (`adjoint(uModelView * transform)`) per instance in the vertex
   shader, a `mat3` varying (WGSL: three flat `vec3`s; 10 varying rows of the 15 / 16 guaranteed): on WebGL2 per
   fragment it cost as if computed per lit sample (the compiler recomputes it at its use in the loop, presumably),
   lit volumes 5-25 % faster (dense blobs 0.86 -> 0.72 ms, the Gaussian volume 1.26 -> 0.95, map-256 0.34 -> 0.33);
   WebGPU level. Measured and dropped: `inverse4` per instance (no gain alone: it runs once before the loop), front
   faces culled by state instead of the shader's `discard` (no gain on either backend). Not done: `flat` for the
   transform varying on GLSL ES 3.0 (equal at all vertices of an instance; not measured).
9. **Optional: a compute marcher on WebGPU** (see Compute), only if the volume pass still dominates after step 5.
10. **Volumes under the transparent objects** (finding 13). Done: every mode draws the volumes (one layer, as before)
    into the target the transparent objects resolve onto, before they do: before `blendedTransparent`, before the
    objects' WBOIT layer (the full-size WBOIT route: a WBOIT layer of their own, evaluated first), before the DPOIT
    peeling; onto the opaque color without postprocessing, into the transparent color target (which postprocessing
    blends over the opaque color) with it. Illumination the same into its transparent target, in all three modes (its
    WBOIT mode no longer puts the volumes into the objects' accumulation). Without transparent objects nothing changed:
    the volumes are composited over the postprocessed image. No extra pass. Transparent objects in front of a volume
    now cover it; those behind it cover it too (the opposite error to before, judged the rarer one with maps). With
    postprocessing the volumes are part of the transparent layer, so its occlusion (`includeTransparent`) and outline
    tint them as they tint that layer. Visual cases: `transparency/*/volume-behind` (a transparent box between the
    camera and the volume) per mode, at the full resolution, with postprocessing, and with illumination.
    - Considered and not done: splitting each ray at the front-most transparent depth (the far part under the
      transparent objects, the near part over them): exact up to that surface, but the transparent depth pass, a second
      ray march (or one writing two targets) and composite, and the far layer placed in every mode; the cost and the
      complexity were judged too high. Ending the rays at the transparent depth: hides the density behind transparent
      surfaces.
11. **Marking and picking for direct volumes** (finding 12), in two parts.
    - The marking variant (done; GLSL and WGSL, `dRenderVariant_marking`): the color variant's ray march (box,
      bricks, opaque depth), its color work replaced by one front-to-back accumulation of the opacity of the samples
      its pass looks at: unmarked ones in the depth pass (`uMarkingType` 1, the occluders), marked ones in the mask pass
      (`uMarkingType` 2; the per-cell or per-instance marker from `tMarker` or `uMarker`). Where it first reaches
      the `surfaceOpacity` param (0.3 by default: edges where the volume looks solid, not where the rays saturate;
      lower for faint transfer functions; `uSurfaceOpacity`, also the pick surface) the ray ends and writes
      that depth and the marker there through `apply-marking` as an opaque fragment (opacity 1); a ray that never
      reaches it discards. One accumulation, not one of all samples and one of the marked: a marked surface behind an
      unmarked one is then occluded by the depth pass (a ghost edge), as for meshes. The ray's phase is constant
      (0.5), not jittered: neighboring rays reach the surface coherently, so its edges are not frayed by the noise, and
      both backends agree.
    - Wiring (done): `getObjectState` allows `marking` for direct volumes, the WGSL program lists it,
      `MarkingPass.planSample` draws `scene.volumes` after `scene.primitives` in both passes, and
      `Scene.markerAverage` counts the visible volumes. The marking passes bind no opaque depth (`tDepth` is far): rays
      cross the whole box and the depth test against the occluders does the rest. Without frag depth (WebGL1 without
      `EXT_frag_depth`) volumes are not drawn in the marking passes: the surface depth can't be written.
    - Picking (done): the `pick` and `depth` variants march to the same surface (all samples count) and write the
      object, instance and cell there (`getDirectVolumeLoci` maps it to `Volume.Cell.Loci`), and its depth (frag depth
      and the packed pick depth; `depth` only for the pick passes without multiple render targets, WebGL only). Fogged
      surfaces are not picked, as other objects (`check_picking_alpha`). `PickPass` draws `scene.volumes` after
      `scene.primitives`, not without frag depth; `getObjectState` refuses only `emissive` and `tracing` now (the
      volumes group is drawn by the volume, marking and pick passes only). A hovered cell is one cell: alone it rarely
      reaches the mask's surface, so a hover shows the label but hardly a highlight unless a granularity extends it
      (e.g. `structure`: the whole volume).
    - Cost (measured, `map-256-highlight-marking`, 2560x1440, everything highlighted): the volume's mask march
      0.39 ms on WebGL2, 0.42 ms on WebGPU, at the mask's full resolution (the volume's color pass at 0.5: 0.33 ms);
      the marking pass as a whole (mask, overlay, its antialiasing) 1.36 vs 0.69 ms a frame (WebGL2), 1.43 vs 0.78
      (WebGPU). The marking layer is cached between camera changes and marking updates, so this is paid while the
      camera moves. Decided: full resolution, not `volumeResolutionScale` (a march at that scale upsampled into the
      mask would save at most about 0.3 ms here, against softer edges).
    - Visual cases (done): `renderer/marking/direct-volume` (half of a volume's cells selected around a highlighted
      opaque box: edges, fill, the box's edges over the volume) and `renderer/pick/direct-volume` (picks a blob, the
      box beside it and the background through `canvas3d.asyncIdentify`, a visual case's `interact` hook, and marks
      what was picked).

Each step: the benchmark cases before and after in this file, the visual cases on both backends, WGSL and GLSL kept in
step.

## Compute (WebGPU)

Compute kernels exist on WebGPU only (`gl/compute.ts`, `GraphicsContext.compute`); WebGL keeps the fragment marcher
either way, so a compute path is a second implementation to keep in step.

- **Ray marching in a compute shader: not first.** The fragment marcher already runs in screen tiles with the
  texture cache's locality; a compute marcher gains from what fragments cannot do: per-tile ray bounds (the tile's
  box entry/exit and opaque depth range in shared memory, tiles outside the box skipped as a whole), ray compaction
  (lanes of terminated rays refilled, against the divergence of early termination) and a stochastic or adaptive
  sample count per tile. It writes a storage texture and needs a composite pass (which step 5 adds anyway). Expected
  gain after steps 1-5 is modest; worth an experiment only if the volume pass still dominates frames then.
- **Preprocessing in compute: yes, where it removes CPU or enables skipping.**
  - Brick min/max per 8³ cells, a reduction kernel over the grid texture: built for the Gaussian volume and dropped
    (its bricks are nearly all visible, step 4); for large maps it could replace the CPU ranges (a quarter of the grid
    texture's build) if that build moves to the GPU too.
  - Brick visibility, or better a distance field over the bricks (each brick: the distance in bricks to the nearest
    visible one, a few passes of a chessboard distance transform), recomputed on transfer changes: the ray jumps
    several bricks at once instead of one. Cheap on the GPU even at 512³ (64 k bricks); on the CPU for WebGL.
  - The grid texture itself (normalized value plus normal, perf-audit V1): from the raw values uploaded as `r32f`,
    instead of 7 `dataOffset` calls per voxel on the CPU; shortens building, not frames.
  - The pre-integration table (step 6), 256x256 entries: small either way, CPU is fine.
- Order: the distance field if bricks one at a time prove the limit for sparse volumes (measure first); the compute
  marcher as an optional step 9 after step 5's measurements.
