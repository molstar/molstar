# Plan: Levels of detail by size on screen

## Status

- **Current**: done. Step 1 measured: `D` 3 looks as `quality` at 720p for half the opaque pass zoomed; with strides
  growing by 10 (1, 10, 100, 1000) the whole view as well, a third of its opaque pass. Decided: `D` in CSS pixels, the
  coarser chain with the pixel presets. Steps 2 (pixel mode, `lodPixelSize`) and 3 (presets, the mesoscale preset, the
  Mesoscale Explorer, the particle representations' default) done.

## Why

- **Measured** (`merging/postsyn-*`, the PetWorld model `postsynstubby_2` as particles with the mesoscale preset,
  RTX 4070 Ti, 1280x720, opaque pass; details in [merged-render-objects.md](merged-render-objects.md), GPU occlusion
  step 1):

  | case                              | WebGPU  | WebGL2  |
  | --------------------------------- | ------- | ------- |
  | whole                             | 4.3 ms  | 5.6 ms  |
  | zoomed (0.3 of the distance)      | 20.4 ms | 25.0 ms |
  | zoomed, LOD preset `performance`  | 5.3 ms  | 6.5 ms  |
  | zoomed, LOD preset `ultra`        | 130 ms  | 159 ms  |
  | zoomed, at 640x360 / 2560x1440    | 15.5 / 29.6 ms | 21.6 / 35.3 ms |
  | zoomed, `ignoreLight`             | 18.3 ms | 22.6 ms |
  | zoomed, `approximate`             | 23.2 ms | 26.7 ms |

  WebGL2 culls on the CPU (no compute), 10k to 13k draws zoomed and 11–12 ms CPU submit against 0.3 ms on WebGPU.
- **The cost is the spheres drawn**: about 0.45 ns a sphere, whole (11–12 M) and zoomed (44 M: levels 0 to 3 4.1,
  16.8, 20.3 and 2.9 M), about 48 spheres a pixel zoomed. Shading is about 2 ms of the 20, the vertex stage alone
  13.4 ms (`sphere-impostors.md` has the per-sphere costs).
- **The levels are chosen by distance in Å** (`Spheres.LodLevelsPresets`: level `j` from `minDistance` to
  `maxDistance`), the same at any resolution, field of view and zoom. What a sphere costs and how it looks depends on
  its size on screen: at 720p and 45° a pixel covers `d * 0.00115` Å at distance `d`, so `quality`'s level 0 ends where
  an atom (radius 1.7 Å) is 3 px, its levels 1 and 2 where their spheres (scaled by `stride^(1 / scaleBias)`) are 1.6
  and 1.3 px; at 2560x1440 they are twice that, in an orthographic view the same for every depth.

## Design

- **A sphere's size on screen**: at view depth `d` a pixel covers `p(d) = d * 2 tan(fov / 2) / H` Å (perspective),
  `p = viewHeight / H` everywhere (orthographic), from the projection (`2 / P[1][1]` is `2 tan(fov / 2)` or the view
  height). `H` is the viewport's height in CSS pixels (divided by the context's `pixelRatio`): a high-DPI display or
  a pixel scale draws the spheres a standard display does, a larger image (`ImagePass`) more. Level `j` (stride `s_j`, scale
  bias `b_j`) draws spheres of radius `r * s_j^(1 / b_j)` (`getLodLevelsValue`'s scale), `r` the sphere radius.
- **Levels by pixel size**: level `j` is drawn where its spheres are at least `D` pixels across, the finest that is:
  where `p(d) < 2 r s_j^(1 / b_j) / D`, and not where the finer level `j - 1` is. `D` (the smallest sphere drawn, in
  pixels) is the one quality knob; the level boundaries follow from the strides and scale biases. Why: one number with
  a visible meaning, the same look at every resolution, zoom and field of view; the presets are strides and scale
  biases only.
- **`r`**: the representation's nominal sphere radius at size factor 1 (`lodRadius`, hidden): 1.7 Å for atoms, a
  carbon's van der Waals radius, whatever their size theme (the presets' `D` calibrated with it); for volume dots 0,
  the largest of their size data (their uniform dot size), so dots of any size switch where they are `D` pixels across
  (`getLodRadius`, again where the sizes changed: the visuals update the sizes before the values). The size factor
  cancels: the stride adjusted for it draws spheres the size of the preset's at size factor 1. A member's own scale (`aSegmentLod`, the stride capped at its sphere count) changes its spheres'
  radius, not its boundaries. Why: the levels stay one list for a merged object (`uLevels`, `uLod` per level) and the
  culling tests one range per cell.
- **Into the existing machinery as distances**: the boundaries in Å per pixel, divided by `k = 2 tan(fov / 2) / H`
  once per frame where `k` changed, are the distance ranges the three tests already use: the cells on the CPU
  (`cullSegment`), the cells on the GPU (`cell-cull`, `uLevels`), the spheres in the vertex stage and the dithered fade
  (`uLod` per level). Orthographic: one level for all, its range `[0, ∞)`, the others empty. Why: no new test, no new
  uniform per sphere; the cost is a few numbers per frame, and the draw lists' keys change only with `k`.
- **`k` from the canvas**, not the pass: the pick pass draws at a lower resolution and must draw the spheres the color
  pass draws; XR views take their own projection's.
- **The fade** between levels (`overlap`) in Å per pixel, converted as the boundaries: adjacent levels fade over the
  same range (the presets have none).
- **At most 4 levels** (`MaxCullLevels`, `aSegmentLod` a vec4), the coarsest open-ended. With a scale bias of 3 a level
  of stride 8 times the finer one covers a doubling of the pixel size; four cover 8x, beyond it the coarsest draws
  spheres smaller than `D` (they shrink with the distance, as now). Measured (step 1): strides growing by 10 (scale
  biases 1, 3, 2.7, 2.3) reach the whole view of the PetWorld model at `D` 3; no more levels needed for it.
- **Params**: `lodPixelSize` (`D`, CSS pixels), 0 for the levels' Å distances (v5 states and apps unchanged); above 0
  the levels' distances are ignored. The level value carries it (`LodLevelsValue.pixel`, boundaries in Å per pixel),
  so merging (its key), the culling and the shaders see one value. The presets get pixel variants (strides and scale
  biases, a `D` each); the mesoscale preset and the Mesoscale Explorer's graphics modes use them.
- **Particles by size on screen**: the particle representations (`ParticleTargetCommonParams`) default to
  `LodPixelPresets.quality` (`lodPixelSize` 3), although particles shipped in v5.11 (a saved state without
  `lodPixelSize` changes look a little, the `quality` look kept). The sphere geometry's default stays 0.

## Steps

1. **Measure the potential** without building it: in the benchmark, distance presets computed from the pixel rule for
   the postsyn camera (`fov`, `H`) at `D` 1.5, 2, 3 and 4; the opaque pass and images against `quality`, zoomed and
   whole, at 1280x720 and 2560x1440. Done (`merging/postsyn-lod-*`, `pixelLodLevels` in the benchmark, its
   `--images`; `quality`'s strides and scale biases, `r` 1.7 Å, 45°; WebGPU, opaque pass):

   | view                | `quality` | `D` 1.5 | `D` 2   | `D` 3   | `D` 4   |
   | ------------------- | --------- | ------- | ------- | ------- | ------- |
   | whole, 1280x720     | 4.28 ms   | 4.29 ms | 4.29 ms | 4.29 ms | 4.29 ms |
   | whole, 2560x1440    | 5.95 ms   | 6.03 ms | 5.97 ms | 5.96 ms | 5.94 ms |
   | zoomed, 1280x720    | 20.5 ms   | 28.9 ms | 18.3 ms | 10.5 ms | 8.1 ms  |
   | zoomed, 2560x1440   | 29.6 ms   | 120 ms  | 78.0 ms | 39.8 ms | 26.1 ms |

   Images against `quality` (pixels differing by more than 16 of 255): zoomed at 720p `D` 3 4.1% (no difference seen),
   `D` 2 18% and 1.5 27% (finer), `D` 4 27% (coarser, blobs visible); the whole view the same at every `D` (all of it at
   the coarsest level, beyond the last boundary). At 2560x1440 `quality` is coarser per pixel (its distances fixed), the
   pixel rule keeps the 720p look and pays for it.

   Coarser chains at `D` 3, 720p (`merging/postsyn-lod-chain-*`), against `quality` (4x4 pixel blocks differing by more
   than 16 of 255, of the model's):

   | strides, scale biases                | whole           | zoomed         |
   | ------------------------------------ | --------------- | -------------- |
   | `quality` (1, 10, 50, 200)           | 4.29 ms         | 20.5 ms        |
   | 1, 10, 50, 200 at `D` 3              | 4.28 ms, 0%     | 10.6 ms, 0.8%  |
   | 1, 10, 100, 1000; 1, 3, 2.7, 2.3     | 1.37 ms, 15.5%  | 7.15 ms, 1.3%  |
   | 1, 10, 100, 1000; 1, 2.5, 2.2, 2     | 1.62 ms, 26.7%  | 11.6 ms, 2.4%  |
   | 1, 20, 400, 8000; 1, 2.5, 2.2, 2     | 0.79 ms, 62.9%  | 8.76 ms, 8.5%  |

   The first coarser chain looks as `quality` (the whole view's speckle a little coarser), the last shows molecules as
   blobs.
2. **Pixel mode** (done): `lodPixelSize` (`Spheres.LodPixelSizeParam`, also the particles' common params); the level
   value in Å per pixel with `pixelSize` (`getLodLevelsValue`, `getLodPixelBoundary`; the preset's distances kept in the
   tuple for the props); `LodView` from the renderer's `update` (the projection, the camera's viewport height over
   `pixelRatio`) to `cull`, `uncull`, `cullSimple` and `prepareCull`; `SegmentedMdbList.prepare` converts the ranges
   (`LodView.getRange`) into the `uLod` uniforms and the ranges the cell tests use, again where the view changed; the
   GPU culler converts per dispatch (`uLevels`, part of what a still frame compares). Tested by `lod-levels`. Measured
   (`merging/postsyn-lod-pixel3*`, the chain of strides 1, 10, 100, 1000 at 3 px, 720p, opaque pass, WebGPU / WebGL2):
   whole 1.36 / 2.04 ms, zoomed 7.12 / 8.47 ms (`quality` 4.28 / 5.6 and 20.5 / 25.0), the images identical to the
   distances computed in the benchmark; orthographic one level for the view (1.33 and 9.04 ms on WebGPU). At 2560x1440
   with a pixel ratio of 1 (twice the CSS pixels) zoomed 35.8 ms.
3. **Presets and apps** (done): `Spheres.LodPixelPresets` per graphics mode, the chain of strides 1, 10, 100, 1000 at
   `D` 6 (`performance`), 4 (`balanced`), 3 (`quality`) and 1.5 (`ultra`), their distances those of the boundaries at
   1280x720 and 45° (for `lodPixelSize` 0); the mesoscale preset (`MesoscalePreset`) and the Mesoscale Explorer's
   graphics modes, presets and LOD controls on them (`lodPixelSize` in `LodParams`). `D` from the chain at 720p, zoomed
   (`merging/postsyn-lod-pixel*`, WebGPU, opaque pass): 1.5 25.9 ms, 2 15.3, 3 7.2, 4 4.7, 6 3.6 (molecules as visible
   blobs); the whole view 1.37 ms at every `D` (all of it at the coarsest level).
