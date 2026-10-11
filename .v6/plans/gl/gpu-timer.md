# Plan: GPU timer with passes as the only timed unit

## Status

- Planning; replaces the earlier WebGPU timer plan (in-pass timestamps, marks charged by their draws). From that plan, steps 0 and 1 are done:
  - the mark audit (table below);
  - the shared timer core `createGraphicsTimer` (`gl/timer.ts`) with a WebGL time-query source;
  - pass boundaries from the contexts, and the `label` option of `beginRenderPass`.
- Step 1 done: pass-only core (`createGraphicsTimer`, one segment per pass, groups charged with the passes that begin inside them, pass leaves, timing of a pass ends with its innermost group), WebGL source with one query per pass, debug placement check; the group calls are still named `mark`/`markEnd` until step 2. Placement warnings in a WebGL Viewer timing run (input for step 2): `Renderer.renderOpaque`, `Renderer.renderBlendedTransparent`, `Renderer.renderPick`, `Renderer.renderMarkingDepth`, `Renderer.renderMarkingMask` and `WboitPass.render` (not in the earlier audit).
- Step 2 done: `group`/`groupEnd` replace `mark`/`markEnd`; the 16 `Renderer.*` groups and `BackgroundPass.render`, `WboitPass.render`, `DpoitPass.render`, `DpoitPass.renderBlendBack`, `TracePass.renderTrace` removed; single-pass groups turned into pass labels (`ColorAccumulate.render` too, instead of moving it); all 109 `beginRenderPass` calls labelled. The debug placement check is silent over the WebGL2, WebGL1 and WebGL2-minimum visual test cases and the Viewer (default, illumination, blended, dpoit). 
- Step 3 done: `timestamp-query` requested where available; `createWebGPUTimerSource` (`gl/webgpu/timer.ts`) adds `timestampWrites` to the pass encoders of timed passes, sums the encoders of a pass across submissions, resolves the used queries at submit into pooled readback buffers (query set and resolve buffer created on the first timed pass); WebGPU render items count draws in `stats.calls` (placement check, `captureStats`). Checked in the Viewer on WebGPU in timing and debug mode: GPU times for groups and passes, no placement warnings; WebGPU visual suite unchanged (33 pass). 
- Step 4 done, over the 33 visual test cases that pass on WebGPU (timing and debug mode, each backend as its own reference): the timing trees have the same shapes (15 of 16 distinct shapes identical; the other is the WebGL-only `hi-Z` group), the same pass labels with the same counts, no placement warnings and no backend errors. GPU times: in the small test renders WebGL reports about 0.05 to 0.1 ms per pass more (time query overhead of ANGLE/D3D11), at Viewer size (1200 × 800) the passes are comparable (e.g. `ssao` 0.41 vs 0.64 ms, `ssao blur` 0.37 vs 0.16 ms, `opaque` 0.27 vs 0.11 ms, `copy` 0.06 vs 0.07 ms on WebGL vs WebGPU). With timing off, no query sets are created and no `timestampWrites` or resolves are encoded. Remaining: step 5 (native compute passes) when they exist; optional quantization note.
- This plan simplifies that core. Related: [webgpu-viewer.md](webgpu-viewer.md), [backend-split.md](backend-split.md), [core-gpu-compute.md](core-gpu-compute.md).

## Goal

One timing model for WebGL and WebGPU:

- **Passes are the only timed unit.** Each render pass (and later compute pass) gets a GPU time from timestamps or time queries at its begin and end. There is nothing inside passes: no in-pass timestamps, no queries switched within a pass.
- **Groups instead of manual marks.** Explicit groups lie clearly outside of passes and sum the passes that began inside them. They don't create GPU work boundaries. Their CPU time is the wall time between begin and end, as today.
- **No special cases.** Tile-based GPUs can't time within a pass by architecture, so no mode exists only for immediate-mode GPUs: WebGL time queries and the Chrome in-pass timestamp feature are not used within passes.
- **No cost.** Timing never splits or adds passes. Only in exceptional cases is a pass split, where the split is logical and cheap; each such split is decided individually and listed below.

## Model

- **Pass**: from `beginRenderPass` to the next `beginRenderPass` or the end of the pass (copies, readback, attachment changes, flush, submit).
  - Its label is the `label` option of `beginRenderPass`, or derived from the target.
  - Timed only while a group is open; charged to the groups open when it began.
  - WebGL: one `TIME_ELAPSED` query per pass. WebGPU: `timestampWrites` on the pass encoders, summed when a pass has several encoders (copies end an encoder, the next draw continues the pass).
- **Group**: `timer.group(label, options)` and `timer.groupEnd(label)`, replacing `mark`/`markEnd` (no aliases).
  - Nested groups form the result tree.
  - GPU time is the sum of the group's passes, or `n/a` when it has none (e.g. CPU-only work and readback).
  - Options as today: `captureStats` and `note`.
  - Its passes are always listed as leaves (`pass <label>`), after its nested groups.
- **Placement rule**: a group begins and ends outside of passes. The pass it times begins inside it (its `beginRenderPass` call is within the group), and no draws of that pass happen after the group ended.
  - In debug mode the timer checks this with the draw counter (`stats.calls`) at group and pass events, and warns when a pass draws while a group is open that was not open at the pass begin, or after a group it was charged to has ended.
  - The rule is a convention checked in debug mode. Groups don't end passes, so timing mode doesn't change the encoded commands.
- **Gaps** between passes (copies, `writeBuffer`, readback) are not timed on either backend.

## Marks to groups (step 2)

After this step there are no marks: the timer API has `group`/`groupEnd` only, every group lies outside of passes, and nothing is timed or labelled inside a pass. Detail below a group comes from labelled passes.

Classification of all 85 `timer.mark` call sites. Sources: the source regions between `mark` and `markEnd` (does the region begin its passes before it draws?), and the debug placement check over the WebGL2 and WebGL1 visual test cases and the Viewer configurations.

### Removed: marks inside a pass

Each of these draws into a pass its caller began. The mark goes away; the caller's `beginRenderPass` gets a `label` naming the work.

| Mark | Calls | Label of the caller's pass (examples) |
| --- | --- | --- |
| 16 `Renderer.*` marks (`renderPick`, `renderDepth`, `renderDepthOpaque`, `renderDepthOpaqueBack`, `renderDepthTransparent`, `renderMarkingDepth`, `renderMarkingMask`, `renderEmissiveOpaque`, `renderEmissiveTransparent`, `renderTracing`, `renderOpaque`, `renderBlendedTransparent`, `renderVolume`, `renderWboitTransparent`, `renderDpoitTransparent`, `renderSolidInteriorCaps`) | the renderer is always called inside a pass | `opaque`, `transparent`, `depth`, `transparent depth`, `opaque back depth`, `marking depth`, `marking mask`, `emissive`, `tracing input`, `volume`, `wboit transparent`, `dpoit transparent`, `pick` |
| `BackgroundPass.render` | one draw in the postprocessing composite pass | the composite pass is labelled `postprocessing` |
| `WboitPass.render` | the evaluate draw in the caller's pass | `wboit evaluate` |
| `DpoitPass.render`, `DpoitPass.renderBlendBack` | one draw each in the caller's pass | `dpoit evaluate`, `dpoit blend back` |
| `TracePass.renderTrace` | the trace draw in an open pass | `trace` |

Solid interior caps and the pick draws of the scene and helpers are part of the pass they share (`opaque`, `pick`, ...); they don't need their own entry.

### Moved: begin just after their pass begins

| Group | Change |
| --- | --- |
| `ColorAccumulate.render` | `group` moves above `webgl.beginRenderPass(framebuffer)` |

Further cases that the debug check finds in step 2 are moved the same way when the move is local, or removed with a label otherwise.

### Turned into labels: groups of exactly one pass

A group around a single pass prints as the group plus one leaf with the same time. These become the label of that pass:

- `CasPass.render` → `cas`; `DofPass.render` → `dof`; `FxaaPass.render` → `fxaa`; `OutlinePass.render` → `outlines`; `ShadowPass.render` → `shadows`;
- `SSAO.opaque` → `ssao`; `SSAO.transparent ` → `ssao transparent` (the trailing spaces of the current labels go);
- `ColorNormalize.render` → `color normalize`; `Grid3dCompute.render` → `grid3d`.

### Kept as groups

They lie outside of passes, and group several passes, calls into code that begins passes, or CPU work and readback:

- frame and helpers: `Canvas3D.render`, `ImagePass.render`, `DrawPass.render`, `StereoCamera.left`/`right`, `MultiSamplePass.renderMultiSample`/`renderTemporalMultiSample`, `PickHelper.render`/`identify`/`asyncIdentify`, `RayHelper.render`/`identify`/`asyncIdentify`, `hi-Z`;
- passes with several passes: `PostprocessingPass.render`, `BloomPass.render`, `SmaaPass.render`, `SSAO.render`, `SSAO.downsample`/`half`/`quarter`/`blurOpaque`/`blurTransparent`, `MarkingPass.redraw`/`updateLayer`/`renderSample`, `DpoitPass.layer`, `IlluminationPass.render`/`renderInput`/`renderMultiSample`, `TracePass.render`/`renderInput`;
- GPU compute and CPU/readback (GPU time `n/a` where no pass begins inside): `GaussianDensityTexture2d`/`3d`, `calcTextureMeshColorSmoothing`, `ColorAccumulate.render` (moved), `Grid3dCompute.renderCumulative`/`renderBatch`/`readPixels`, `createHistogramPyramid`, `getHistopyramidSum`, `calcActiveVoxels`, `createIsosurfaceBuffers`, `extractIsosurface`, `createGaussianSurfaceTextureMesh`, `createStructureGaussianSurfaceTextureMesh`, `PickBuffers.read`/`asyncRead`, alpha-orbitals `createSphericalCollocationGrid`/`createSphericalCollocationDensityGrid`.

### API

- `GraphicsTimer.group(label, options)` and `groupEnd(label)` replace `mark`/`markEnd`; `TimerOptions` stays (`captureStats`, `note`). No aliases; `mark`/`markEnd` are removed, recorded in [breaking-v6-changes.md](../breaking-v6-changes.md).
- Callers: the graphics package, the plugin (none today besides `resolve`), the alpha-orbitals extension.
- The WebGPU context's stub timer gets the same interface (`createGraphicsTimer(undefined, stats)` until step 3).

### Done when

- No `mark`/`markEnd` in the code base, no group inside a pass: the debug placement check is silent over all visual test cases (WebGL2, WebGL1, WebGL2 minimum) and the Viewer configurations (default, illumination, blended, wboit, dpoit, bloom, stereo).
- Timing output: the group trees of a Viewer timing run are the old trees without the removed and relabelled marks, with labelled pass leaves; no `n/a` except for groups without passes.
- Visual profiles byte-identical (labels and groups don't affect rendering).

### Splitting passes (exceptional)

Considered and rejected, because each split adds a load/store of the attachments:

- opaque and blended transparent draws in one color pass;
- emissive opaque and transparent draws in one pass;
- the background with the postprocessing composite.

None is planned. A split is only considered when the two parts are logically separate passes anyway and the extra load/store is small (e.g. a small target), and is then recorded here.

## Core simplification

Relative to the current `createGraphicsTimer`:

- **Remove** in-pass sources:
  - `TimerSource.inPassTimestamps` and the segment switches at marks;
  - charging by draws (deepest mark containing the draws), `insidePasses` in results, and the "inside pass" `NaN` path.
- **Keep**: pass records and labels, segment per pass, groups charged with the passes that began inside them, pass leaves (always shown, the `passes` option is removed), stats capture, moving averages, dropping invalid (disjoint) results, waiting for passes that end after their group.
- **Add**: the debug-mode placement check (draw counter at group and pass events), and `group`/`groupEnd` names.
  - `mark`/`markEnd` are removed; callers in the graphics package, the plugin and the extensions (alpha-orbitals) switch to the new names; recorded in [breaking-v6-changes.md](../breaking-v6-changes.md).
- **WebGL source**: one query per pass, begun at the pass begin and ended at the pass end (sequential, never nested), no queries outside passes.

## Steps

1. **Core**: simplify `createGraphicsTimer` as above; WebGL source with one query per pass; debug placement check; unit tests (charging by pass begin, nested groups, passes ending after their group, placement warnings, stats capture, invalid results).
2. **Marks to groups** (details in "Marks to groups (step 2)"):
   - rename the timer calls;
   - remove the `Renderer.*`, `BackgroundPass.render` and `TracePass.renderTrace` marks;
   - add pass labels in the callers;
   - fix groups that violate the placement rule (the debug check over all visual test cases and the Viewer configurations must be silent).
3. **WebGPU source**: `timestamp-query` in `OptionalFeatures`; `timestampWrites` on pass encoders while a group is open; query set and readback buffer pools; encoders summed per pass.
4. **Check**: timing runs on both backends print the same group trees with the same pass labels; per-pass GPU times plausible against each other; no validation errors; frames with timing off encode the same commands as before.
5. **Native compute passes** (when added): timed like render passes, through `timestampWrites` on compute pass descriptors.

## Validation

- WebGL and WebGPU visual profiles unchanged (timing is off in visual tests; labels don't affect rendering).
- The debug placement check is silent over all visual test cases and the Viewer configurations (default, illumination, blended, wboit, dpoit, bloom, stereo).
- Timing output: the group trees are the old trees without the removed marks; the removed details appear as labelled pass leaves.

## Decisions

- `mark`/`markEnd` are renamed to `group`/`groupEnd`, without aliases.
- Pass leaves are shown by default; there is no `passes` option.
