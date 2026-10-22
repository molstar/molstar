# Animation playback: per-frame cost

Goal: make the per-frame work of playback proportional to what the frame changes. Two kinds of playback are covered:
state-snapshot and MVS animations, and trajectory (model-index) playback. The load side (building MVS frames, the
state tree, representations) is covered by [mvs-optimizations](../mvs-optimizations/research/README.md) and
[perf-audit](../perf-audit/README.md); this is what happens once playback runs. Branch: `v6-animation-playback`.

The findings below come from two read-only traces of master at `93e710836` (2026-10-22). They are claims until the
bench (step 0) or a test confirms them; ids are `S` for snapshot/MVS playback and `T` for trajectory playback. Paths are
relative to `packages/`.

## Status

| Step | Ids                    | What                                                                      | Status |
| ---- | ---------------------- | ------------------------------------------------------------------------- | ------ |
| 0    | —                      | Playback bench: ms and bytes allocated per frame, both kinds              | done: [Baseline](#baseline-step-0) |
| 1    | S1–S5, S9, B2          | Snapshot playback: the per-tick re-apply bug and per-tick waste           | done: S1–S4, S9, B2; S5 dropped |
| 2    | T1, T2, B1             | Trajectory: the model cell and the shared root (correctness first)        | done: T1, B1 (with T5's focus fix); T2 with T3 |
| 3    | T3, T4                 | Trajectory: remap instead of recreate; keep the topology caches           | done: T3 (and T2); T4 dropped |
| 4    | R1                     | Representations: no theme or params rebuild when nothing they read moved  | done: the themes; the param definition cache removed (review gate 2); multilayer selectors open |
| 5    | T5                     | Behaviors: selection, focus, the empty update, Prism marks                | done: the empty update, SelectLoci, the focus (with B1); Prism marks kept |
| 6    | T6                     | Visuals: unchanged buffers stay, complex visuals keep colors              | dropped (measured; open for a browser) |
| 7    | S6, S7                 | A frame applied as a diff of the changed refs; unchanged children skipped | done: S7 in the effect transforms, bundle layers; S6 dropped (2.3% of a frame) |
| 8    | S8                     | State interpolation: equality and the ref match                           | done |

Order: the bench, then the bugs (1, 2), then the work per frame from the bottom of the tree to the top (3–6), then the
state-level change (7). Step 7 changes how a snapshot is applied, so it needs the core state spec, and the steps before
it make its gain measurable on its own.

## Summary: baseline → the branch

The bench at the end (median of 3 runs, the MVS and snapshot scenes in one process and the trajectory scenes in
another), against the baseline below. ms/applied: per frame that changed the state; ms/tick: the median tick (ticks
between frames included).

| Scene                          |    ms/applied | ms/tick med |        KB/tick | Cells/tick u / R            | setSnapshot |       Log |     Marks |
| ------------------------------ | ------------: | ----------: | -------------: | --------------------------- | ----------: | --------: | --------: |
| snapshots:plain                |     1.1 → 0.9 | 0.03 → 0.00 |       98 → 47  | 0.05 / 0 → 0.05 / 0         |   181 → 31  |    9 → 0  |         0 |
| mvs:exosome                    |     4.7 → 3.0 | 3.98 → 2.78 | 3,946 → 2,724  | 17.52 / 0 → 17.58 / 0       |  601 → 419  | 10,515 → 0 |         0 |
| traj:nmr:preset                |     2.2 → 1.8 | 0.45 → 0.03 | 1,347 → 926    | 2.25 / 2.25 → 4.5 / 0       |           0 |         0 |   27 → 0  |
| traj:nmr:preset+focus+select   |     5.8 → 4.7 | 2.63 → 0.03 | 6,435 → 1,985  | 5.7 / 2.65 → 6.75 / 0       |           0 |         0 |       108 |
| traj:synth:preset              |    13.2 → 9.2 | 0.16 → 0.03 | 8,919 → 4,068  | 2.44 / 1.95 → 4.39 / 0      |           0 |         0 |  117 → 0  |
| traj:synth:preset+focus+select |   14.5 → 12.7 | 9.07 → 0.03 | 10,746 → 5,303 | 6.3 / 2.42 → 6.83 / 0       |           0 |         0 |       468 |
| traj:synth:cartoon             |    11.9 → 8.4 | 0.15 → 0.03 | 6,689 → 3,629  | 1.46 / 1.46 → 2.92 / 0      |           0 |         0 |   39 → 0  |
| traj:synth:ball-and-stick      |     8.0 → 5.1 | 0.17 → 0.03 | 6,266 → 1,361  | 1.46 / 1.46 → 2.92 / 0      |           0 |         0 |   39 → 0  |
| traj:synth:gaussian-surface    | 296.0 → 292.6 | 0.23 → 0.03 | 11,699 → 11,381 | 1.46 / 1.46 → 2.92 / 0     |           0 |         0 |   39 → 0  |

New scenes (no baseline): snapshots:effects 0.72 ms/applied; traj:nmr:preset+focus 4.64; traj:synth:preset+overpaint
9.17. With a selection the marks stay (they are applied again on every frame: review gate 1). What is left in a frame:
the geometry for the new coordinates (trajectories), and on the exosome the `gl/` checks and merging (Open).

## Baseline (step 0)

`tests/src/pipeline/playback.ts`, run by `node scripts/pipeline-bench.mjs playback` (its header has the options). A real
`PluginContext` (the default spec and the MVS behavior) with the stub canvas of `many.ts`, now shared
(`canvas-stub.ts`), driven by `managers.animation.tick(t, true)` on a synthetic 60 Hz clock: every tick awaits its
frame, so a tick costs what the browser's loop pays when it keeps up. Trajectories play `AnimateModelIndex` at 30
frames/s (`computed`), so every other tick has no new frame; snapshots play `AnimateStateSnapshotTransition`.

- **Columns.** Applied: ticks with any cell event. ms: per tick, median and max, and the median of the applied ticks.
  KB/tick: V8 `total_allocated_bytes` over the range, divided by the ticks; one read per range, because a read itself
  allocates ~165 KB. u / R: cells updated in place and recreated (`object.updated`), per tick. Camera: resets, and
  restarts (a reset with a duration while the previous transition still runs). setSnapshot: whole-tree
  `State.setSnapshot` calls. `changed`: the snapshot manager's `changed` emits. Marks: `canvas3d.mark` calls.
- **Scenes.** snapshots:plain is 1cbs as a transition of 30 frames of 100 ms over three looks (preset, chain-id,
  element colors), with a 400 ms camera move. mvs:exosome is snapshot 2 of the exosome story (1,022 frames, 27 s; the
  first 600 ticks), from `mol-view-stories/.../examples/exosome/story.mvsx` (local only: `--mvsx <file>`; motm1 needs
  the PDBe files, so it is not in the table). traj:nmr is 1grm (`data/examples`, 5 models of 272 atoms, the mmCIF
  path). traj:synth is 7qpd's topology with 20 frames of jittered coordinates (11.5k atoms, the topology + coordinates
  path, as a DCD). Representations: the default preset (components, cartoon and ball-and-stick), or one representation
  on everything; +focus+select adds a focused residue and a selected chain.
- **Machine.** Apple M4, Node 25.6, master `93e710836` plus the bench; median of 3 runs.

| Scene                          | Ticks | Applied | ms/tick med (max) | ms/applied med | KB/tick | Cells/tick u / R | Camera (restarts) |    Log | setSnapshot | `changed` | Marks |
| ------------------------------ | ----: | ------: | ----------------: | -------------: | ------: | ---------------: | ----------------: | -----: | ----------: | --------: | ----: |
| snapshots:plain                |   181 |       3 |         0.03 (2) |            1.1 |      98 |         0.05 / 0 |             7 (6) |      9 |         181 |       179 |     0 |
| mvs:exosome                    |   600 |     417 |        3.98 (18) |            4.7 |   3,946 |        17.52 / 0 |             0 (0) | 10,515 |         601 |       599 |     0 |
| traj:nmr:preset                |    20 |       9 |         0.45 (4) |            2.2 |   1,347 |      2.25 / 2.25 |                 0 |      0 |           0 |         0 |    27 |
| traj:nmr:preset+focus+select   |    20 |      17 |        2.63 (18) |            5.8 |   6,435 |       5.7 / 2.65 |                 1 |      0 |           0 |         0 |   108 |
| traj:synth:preset              |    80 |      39 |        0.16 (34) |           13.2 |   8,919 |      2.44 / 1.95 |                 0 |      0 |           0 |         0 |   117 |
| traj:synth:preset+focus+select |    80 |      77 |        9.07 (26) |           14.5 |  10,746 |       6.3 / 2.42 |                 1 |      0 |           0 |         0 |   468 |
| traj:synth:cartoon             |    80 |      39 |        0.15 (23) |           11.9 |   6,689 |      1.46 / 1.46 |                 0 |      0 |           0 |         0 |    39 |
| traj:synth:ball-and-stick      |    80 |      39 |        0.17 (23) |            8.0 |   6,266 |      1.46 / 1.46 |                 0 |      0 |           0 |         0 |    39 |
| traj:synth:gaussian-surface    |    80 |      39 |       0.23 (333) |          296.0 |  11,699 |      1.46 / 1.46 |                 0 |      0 |           0 |         0 |    39 |

What the baseline confirms and what it adds:

- **S1, S2, S3 hold.** Every tick runs a whole-tree `setSnapshot` (181 of 181, 601 of 600), and the manager emits
  `changed` on each. During frame 0 (100 ms, 6 ticks) the camera is reset 7 times and the transition restarted 6 times.
  The exosome logs 10,515 entries in 10 s.
- **The trajectory recreates cells** on every frame (R, about 3 per frame with one representation, 4 with the preset):
  T1 and T3. With a focus and a selection, the ticks without a new frame also update (77 of 80 ticks applied, not 39):
  T5.
- **Profile of the exosome frames** (`--profile`, 2.1 s sampled over 600 ticks): the tree diff and `fromJSON` are not
  where the time goes at this size (75 nodes). Shares of the sampled time:
  - render-object merging in `prepareUpdate` (`repr/merged.ts` over `gl/render-object.ts`): 22%;
  - `relinkCells` → `getEffectiveDependsOn` → `tryGetDefaultParamDefinition`: 14%. For each scheduled cell, it builds
    the transformer's whole param definition (the representation's `getApplicableTypes`) only to find its refs. This
    one is new: **S9** (below);
  - the representation's `params` in `resolveParams`: 10%, and theme creation: 7% (R1).
- **New findings:**
  - **B1.** With a focus, playing 1grm throws in the focus representation's interactions
    (`interactions-intra-unit-cylinder.ts:66`, `InteractionsProvider.get(structure).value` is undefined after a frame).
    To be checked with T2 and T5.
  - **T5, more.** The trajectory animation updates each model cell with a recipe that returns a new object
    (`trajectory.ts:169`), so every tick goes through `produce`, and mutative warns "The return value does not contain
    any draft".

## Results

One line per id as it lands: the commit's before → after on the bench (median of 3 runs, same machine as the baseline).

- **S1** (state-snapshots.ts: the skip reads `currentTransitionFrame`, which the animations now write). plain: 181 → 30
  `setSnapshot` per 181 ticks, camera resets 7 (6 restarts) → 1 (0), 98 → 43 KB/tick. exosome: 601 → 418
  `setSnapshot` per 600 ticks; ms/tick within noise (3.98 → 4.32 med), since the skipped applies changed no cell. Test:
  `plugin/core/src/state/animation/_test/state-snapshots.test.ts` (36 and 30 applies before, 5 and 4 after).

- **S2** (`State.setSnapshot(snapshot, options)`; the playing and the scrubbed frames pass `doNotLogTiming`). exosome:
  10,515 → 0 log entries per 600 ticks, 4.32 → 3.12 ms/tick med, 5.33 → 4.63 ms/applied; plain: 9 → 0 entries.
  Tests: two in `state-snapshots.test.ts` (12 and 8 entries before).

- **B2** (`getStateTransitionFrameIndex`: frame `i` from the end of the frames before it until its own end). A fix of
  what is shown, not of the cost; S1's test now applies 6 and 5 frames (frame 1 included). Test:
  `plugin/core/src/state/_test/transition-frames.test.ts` (failed before: the hold of a 100/1000/100 ms transition
  showed frame 2). Also in 5.x: in the root `CHANGELOG.md` with S1 and S2.

- **S3** (`events.animationTimeChanged` for the time; `changed` for the entries, the current one, playing). exosome:
  416 → 0 `changed` per 600 ticks; ms within noise (4.23 med). The transition now sets the time before its skip, so
  the time also moves while a frame holds (after S1 it moved on frame changes only). Readers of the time: the classic
  viewport controls and Prism's transition source; the others (the snapshot list, two markdown descriptions,
  mvs-stories' markdown, video export) no longer wake per tick. Breaking: in `.v6/plans/breaking-v6-changes.md`.
  **Dropped: the binary search.** The scan is 0.38 µs per call over 1,022 frames (0.06 µs binary), two calls per tick
  against 3–4 ms. Test: in `state-snapshots.test.ts` (failed before: no event).

- **S4** (`state.ts`: the transformer ids and provider names a snapshot uses, read once per tree and per frames array,
  `usedNames`). The checks of exosome's snapshot 2 (1,022 frames) on each apply: 18.75 → 0.02 ms (the bench's MVS
  scene prints them). Same warnings in the same order. Tests: in `report-unregistered-names.test.ts` (100
  `StateTransformer.has` calls for two repeat checks before, 2 after).

- **S9** (`transform.ts`: the default param schema cached per transformer, not per transform). exosome: 4.62 → 3.72
  ms/applied, 4.23 → 2.94 ms/tick med; plain 0.94 → 0.81 ms/applied. Test:
  `core/src/state/_test/snapshot-frames.test.ts` (4 schema builds for 4 frames before, 1 after).

- **S5: dropped.** It holds in a test: a child whose serialized transform the frames share is normalized again on each
  frame, because `fromJSON` makes a new `Transform`. A per-state memo of the read transforms fixed that. But on the
  exosome it changes nothing: `normalizeParams` ran 10,520 times over 300 ticks with and without it, and ms/applied
  stayed at 3.72 → 3.74. The nodes a frame evaluates there have new versions. Not worth a new rule in `tree.md`.

- **T1** (`ModelFromTrajectory.update`: the new frame's model in place, nothing disposed; `dispose` on delete stays).
  Confirmed first: a frame change disposed the assets of a property attached to any frame. ms/applied: synth preset
  13.17 → 9.53, cartoon 11.94 → 8.35, ball-and-stick 8.01 → 5.05, preset+focus+select 14.47 → 10.85; nmr+focus 5.75 →
  3.40. KB/tick: synth preset 8,919 → 4,351, ball-and-stick 6,266 → 1,443. Recreated cells per frame: one fewer (3 →
  2 with one representation). Gaussian surface unchanged (296 → 301, the surface itself). Tests:
  `plugin/core/src/state/animation/_test/trajectory.test.ts` (the asset disposed, the cell recreated before). Also in
  5.x: root `CHANGELOG.md`. After the review gates: another trajectory (the data read again) recreates the cell (gate
  1); a frame that does not share the properties (another topology: only `sameAsPrevious` frames share
  `{ ...previous }`) disposes the old ones, and `CustomModelProperties` releases its references on the model that held
  them; `ModelWithCoordinates` below builds again from the new model (gate 2).

- **T2: not reached today, done with T3.** A static component never gets to its remap path during playback:
  `Structure.areEquivalent(a, cache.source)` compares unit conformations by identity, so every frame recreates the
  component and runs its query again (T3). Once T3 remaps the component, the separate roots of `remapModel` would
  appear, so T3's commit passes the shared root. The test (`trajectory.test.ts`, "share the remapped root") passes
  before and after.
- **T5 (focus) and B1** (`structure-focus-representation.ts`: `focus()` commits nothing when the target's residue bundle
  is the one its selection has; the param changes force it). Per frame the focus manager remaps its loci and the
  behavior committed a second update, which rebuilt the surroundings expression, recreated that selection and updated
  the three focus representations, after the frame's own update had already moved them. B1 was that second update
  still queued when the plugin is disposed: it rebuilt the interactions on a structure whose properties the dispose
  had cleared. synth+focus+select: 77 → 39 ticks with updates, 7,028 → 5,727 KB/tick, 2.16 → 0.14 ms/tick med (the
  frame's work is now in one tick: 10.85 → 11.84 ms/applied); nmr+focus+select 2,547 → 2,014 KB/tick. Tests: two in
  `trajectory.test.ts` (2 updates per frame, and the throw after dispose, before).

- **T3 and T2** (no cell is recreated on a trajectory frame). `CustomStructureProperties` updates in place on a new
  structure: it disposes the old structure's own properties and attaches the new one's. A static component remaps
  when the new source has the same hierarchy (its query reads topology only, and remapped units keep their bonds) onto
  `a.root`, the root already remapped (T2: `Structure.remapModel(m, parent?)`, model `_spec/structure.md` §3). A
  whole-source component takes the new source. The effect transforms (overpaint, transparency, substance, emissive,
  clipping, wiggle) put the layers on the new structure in place. Recreated cells per tick → 0 in every scene (from
  0.97–2.92). ms/applied: synth preset 9.69 → 9.22, +overpaint 9.74 → 9.20, cartoon 8.75 → 8.29, ball-and-stick 5.08
  → 5.07; nmr+focus+select 6.44 → 4.48. KB/tick: synth preset 4,989 → 4,139, +overpaint 4,772 → 4,161, nmr+focus+select
  5,258 → 1,936. Measured before and after back to back. Tests: `trajectory.test.ts` (3 recreates per frame before;
  without the shared root the T2 test fails once components remap). Changed: a whole-source static component stays its
  source on a new frame (it was a remapped copy), in `structure-component-element-schema.test.ts`. After the review
  gates, the remap also needs equal unit operators (gate 1: a frame that also moves the structure) and no dynamic bonds
  (gate 2: the "plus connected" queries follow bonds computed per frame); otherwise the query runs again.
  **Dropped: expression components** (`structure-component.ts`, `expression`/`script`). An expression may read
  coordinates (`includeSurroundings`, distances), so running it again is the only safe answer without analyzing the
  expression.

- **T4: dropped.** Profiled after T3 (synth preset, 50 frames; and 1tup-flat-32 with ball-and-stick, 288 units,
  186k atoms): `unitSymmetryGroups`, `hashCode`, `models` and `carbohydrates` of the remapped structure are not in the
  top 60 entries. The frame is the new geometry (the cartoon mesh 77%; bond cylinders and spheres 70% with
  ball-and-stick) and the unit boundaries (15–20%), which read coordinates. A unit's `transientCache` (valence model,
  interaction features) reads coordinates too, so nothing there can be carried over.

- **R1 (themes)** (`representation.ts` update: the color and size themes are kept when the representation's own `type`
  params changed and the structure and the themes' params did not, with their dependencies' references). After review
  gate 1, an update with the params as they were (an ancestor changed, e.g. the secondary structure type of the
  structure's properties) makes the themes again: a theme holds the version of a property it reads. exosome: 3.72 → 2.95 ms/applied, 3,632 → 2,812 KB/tick, 2.94 → 2.40
  ms/tick med. Trajectories are unchanged by it: a new structure makes the themes again. Tests:
  `plugin/core/src/state/transforms/_test/representation-update.test.ts` (a theme made on an alpha change before).
  Method note: a trajectory scene measured in the same process after the exosome scene reads up to 40% slower (the
  exosome allocates 2.3 GB first), so trajectory numbers come from runs of their own.

- **R1 (param definition): removed after review gate 2.** It was: the definition cached by structure, checked against a
  new `version` of the representation and theme registries. Stale when a theme's `getParams` reads property data its
  own `ensureCustomProperties` attaches in the update (the PDBe structure quality report's issue types); the cache and
  the `version` are gone, the 6% below with them. Measured while it was in: exosome: 2.95 → 2.76 ms/applied, 2,812 →
  2,411 KB/tick. Trajectories unchanged (a new structure per frame). Tests: in `representation-update.test.ts`
  (`getApplicableTypes` 3 times for 3 updates before; a registration builds it again).
  **Not done:** `computeChildAwareTransformGroups` is in none of the profiles. **Open:** the MVS multilayer color
  theme's per-layer selectors could be cached by structure and selector. It matters for color animations only, and the
  fixtures that have them (features, motm1) need PDBe files, so it can't be measured here.

- **T5 (the empty update)** (`animation/trajectory.ts`: a cell is updated only when its frame index changes, with the
  params object rather than a recipe; nothing is committed when none changes). The ticks within a frame: 0.12 → 0.03
  ms/tick med in every trajectory scene; `isUpdating` toggled 61 → 2 times over 60 ticks of 3 frames; mutative's
  warning per tick is gone. Test: in `trajectory.test.ts` (61 toggles before).

- **T5 (SelectLoci)** (`behavior/dynamic/representation.ts`: no re-marking when the updated structure has the same
  hierarchy and nothing selected; with a selection the marks are applied again, as before). Narrowed after review gate
  1: an operator change resets a units visual's markers (`updateTransform`), and bonds or contacts computed per frame
  (dynamic bonds, interactions) move the groups the marks are on, so with a selection only re-marking is right.
  `canvas3d.mark` calls without a selection: synth preset 117 → 0 over 80 ticks, cartoon 39 → 0, nmr preset 27 → 0,
  nmr+focus 54 → 8; with one, unchanged (synth+focus+select 468). Time within noise in Node, where the stub canvas
  makes a mark cheap; in the browser a mark rewrites the marker data and walks the scene. Tests: in
  `trajectory.test.ts` (2 mark calls per frame without a selection before; with one, marked again).

- **T5 (Prism marks): kept.** Correctness needs it. A frame updates a representation's own overpaint cells in place
  (T3), and they set its overpaint to their layers only. The marks paint over those, so they paint again once per
  animation frame, on the changed representations only. Not measured: Prism, not the Node bench.
- **T6: dropped, open for a browser measurement.** The bench now prints the KB of render-object arrays whose value cell
  changed per applied frame, what the renderer would upload again: cartoon 14,127, ball-and-stick 4,831, gaussian
  surface 31,152 (synth). In the cartoon, indices and groups are 6.8 of the 14.1 MB. They are equal across frames, but
  the builders write the new frame into the old arrays (`ChunkedArray.create` with the cell's array as the first
  chunk), so the cell gets the same array back and must update. A content compare finds nothing (tried: unchanged on
  every scene). Saving the upload needs the builders to allocate fresh topology arrays and compare them: CPU and
  allocation for upload, which only a browser can weigh. `ComplexVisual`'s forced color and size update is in none of
  the profiles.
- **Open, `gl/` (frozen; for the WebGPU work).** Of an exosome frame after step 4: `prepareUpdate` →
  `hasForeignResources` 22% (a walk over every value cell of the render object and its merged members, on every visual
  update, for the backend switch), and render-object merging (`repr/merged.ts` → `canMergeRenderObjects`) 16%. The
  callers are in `repr/`, but what they check is the backend-switch logic.

- **S8** (`state-interpolation.ts`: targets by ref in a map; a node is updated when its interpolated params differ from
  the cell's params or from what the animation set last; nothing is committed on a tick that changes nothing). The
  comparison was of the cell with the params (`shallowEqual(cell, params)`, always false), so every differing
  transform was updated on every tick. Over 58 ticks of two 4-node snapshots: 58 → 1 commits, the same 4 nodes
  changed. Not in the bench: an experimental animation with no scene. Test:
  `plugin/core/src/state/animation/_test/state-interpolation.test.ts` (58 commits before).

- **S7** (in the effect transforms, as decided: overpaint, transparency, substance, emissive, clipping and wiggle skip
  `ofScript`/`ofBundle` in `update` when the layers' params and the structure are those they last ran on; no
  reconciler or spec change). The bench gets `snapshots:effects` (1cbs, a script overpaint and transparency below each
  representation, frames changing the alpha): 1.00 → 0.74 ms/applied, 55 → 38 KB/tick. MVS does not use these
  transforms (its opacity and emissive are representation params), so the exosome is unchanged. Test:
  `plugin/core/src/state/transforms/_test/effects-update.test.ts` (one `ofScript` per update before). **Narrowed after
  review gate 2 to bundle layers:** a script may read a custom property that an ancestor's update computes again on the
  same structure (ASA's `is-buried`), so script layers run again. The scene's overpaint is now a bundle (its
  transparency a script): 0.86 → 0.80 ms/applied, 55 → 42 KB/tick.
- **S6: dropped**, by the rule decided before step 7. `fromJSON` + `findDeletes` + `findUpdateRoots` + `syncNewStates`
  are 39.5 of 1,709 ms sampled over 600 exosome ticks after steps 1–5: 2.3% of a frame, under the 10% line.
  `reconcile.md` is unchanged.

## Review gate 1 (after step 6)

A fresh reviewer read `master...23940ded6`. Findings and what was done:

1. **SelectLoci dropped the marks after a transform** (bug). `areHierarchiesEqual` ignores operators, and an operator
   change resets a units visual's markers. Fixed: no re-marking only without a selection (`5d9f852f3`).
2. **A static component kept stale operators when remapped** (bug): a frame that also moves the structure. Fixed: the
   remap needs equal unit operators (`e3b036401`).
3. **Kept themes went stale when a dependency was computed again** (bug): an ancestor's update (the structure's
   properties) evaluates the representation with equal params, while the theme holds the old property version. Fixed:
   the themes are kept only on the representation's own update (`d109014d0`).
4. **Stale marks where bonds or contacts change per frame** (risk: dynamic bonds, interactions). Covered by the fix
   of 1: with a selection, marking again.
5. **ModelFromTrajectory leaked a new trajectory's predecessor's properties** (risk): a re-read trajectory kept the
   cell in place. Fixed: another trajectory recreates the cell (`c1952c0f6`).
6. **The per-transformer default schema is computed once per process** (nit). Kept: `params(undefined, undefined)`
   has no source and no context, and the reviewer found no transformer whose default schema reads a registry.
7. **S8 skips a cell whose interpolated params equal what it set last**, so a change made by someone else in the
   middle of the animation is not reset (nit). Kept: the animation does not fight an edit.

The reviewer checked as correct: the snapshot skip, the camera, the frame index function, the time event's readers,
`CustomStructureProperties` in place, the effects in place, the shared root, the focus skip, the registry versions,
the snapshot names cache, the trajectory apply loop, and that the other tests and the model `_spec` edits are right.

## Review gate 2 (after step 7)

A fresh reviewer read `master...5fd5a1acb`. Findings and what was done:

1. **S7 left script effects stale** when a script reads a custom property recomputed on the same structure (bug).
   Fixed: the early return is for bundle layers only (`8f5bb4603`).
2. **The R1 param definition cache went stale** when a theme's params read property data (bug: the PDBe quality
   report). Fixed by removing the cache and the registries' `version` (`20c00fd0a`).
3. **Frames that do not share custom properties** (another topology) kept their assets and references (risk). Fixed:
   disposed when not shared, references released on the model that held them (`22ae67ff5`).
4. **Static components with dynamic bonds** kept frame 0's connected elements (risk). Fixed: the query runs again
   (`22ae67ff5`).
5. **`ModelWithCoordinates` below a model updated in place** kept the old model (risk). Fixed: it builds again from
   the new one (`22ae67ff5`).
6. **SelectLoci re-marked on updates with the same data** (nit). Fixed: those keep the marks again (`22ae67ff5`).
7. **Two tests checked less than their names** (nit). The re-marking test now checks that the marks are on the new
   frame's structure; the S7 test covers script layers.
8. **A transition ended without its last frame when that slot is shorter than a tick** (not a regression). Fixed:
   the end applies the last frame unless it is shown (`22ae67ff5`).

The reviewer checked as correct: the snapshot skip and camera, the frame index function, the time event, the S9
cache, the names cache, the trajectory identity check, `CustomStructureProperties`, the static remap and shared root,
the effects in place, the kept themes, the focus skip, the trajectory loop, S8, and the spec and breaking-change edits.

## Snapshot and MVS playback

- **S1, bug (checked in the code).** `plugin/core/src/state/animation/built-in/state-snapshots.ts`: both animations skip
  a tick when `frameIndex === animState.currentTransitionFrame` (:143, :239), but they only ever set
  `currentAnimationFrame` (:157, :165, :258), so the skip never fires. Every animation-loop tick then runs a full
  `setSnapshot`. During frame 0 it also re-applies the camera and canvas props and calls `requestCameraReset`, which
  restarts the camera transition on every tick: the move is stretched, not only wasted.
- **B2, frames shown one slot early (found with S1's test).** `getStateTransitionFrameIndex` returns `i` for the first
  `i` whose start (the sum of the durations before it) is at least `t`, so during frame `i`'s slot it returns `i + 1`.
  Frame 1 is never shown, and during an MVS hold (a long frame) the frame after it is shown. Scrubbing
  (`setSnapshotAnimationFrame`) has the same offset. Fixed with S3, which replaces the scan.
- **S2, log growth (checked).** The snapshot animations do not pass `doNotLogTiming`, unlike the trajectory, spin,
  explode, unwind and interpolation animations. So every frame logs "Updated X in N ms" for each updated cell, without
  limit (`context.ts:692`), and re-renders the log view.
- **S3, per-tick waste.** `getStateTransitionFrameIndex` (`plugin/core/src/state.ts:520-537`) scans the frames linearly,
  twice per tick: use a binary search over cumulative durations. `updateState({ currentAnimationTimeMs })` fires the
  snapshot manager's `changed` on every tick, which wakes the snapshot UI and Prism's snapshots source
  (`prism/viewer/src/features/snapshots/source.ts:149-153`): give the time its own event or throttle it.
- **S4, per-snapshot waste.** `validateSnapshotTransformers` and `reportUnregisteredNames` (`state.ts:168, 206`) walk
  every frame's whole tree on each `setSnapshot`: O(frames × transforms), about 2.9k × 750 on the exosome. Do it once
  per loaded snapshot.
- **S5, params normalized on every frame.** `resolveParams` (`core/src/state/state.ts:1373-1387`) keeps
  `_normalized_param_version` on the `Transform`, which `StateTree.fromJSON` creates fresh every frame, so it never
  matches. Key it by something that survives, e.g. the params object, which frames share by reference.
- **S6, a frame is applied as a whole tree.** Each frame holds the full serialized tree. Applying it means a
  `fromJSON` (new immutable maps, a new `Transform` per node) plus three walks over the whole tree (`findDeletes`,
  `findUpdateRoots`, `syncNewStates`), even though consecutive frames share transform objects since MVS load steps 3–4.
  Compare a frame's transforms with the previous frame's by object identity (version as the fallback), and apply only
  the changed refs as a builder, so that the reconcile's fast track and `editsOnly` apply. When the tree's shape
  differs (adds, deletes, moves), fall back to the full `setSnapshot`. This is the neighbour of MVS research step 5,
  "frames computed at playback": say in the PR how they relate.
- **S7, unchanged children re-run.** `updateSubtree` (`core/src/state/state.ts:1364-1370`) recurses into every
  descendant of an updated cell. For example, an overpaint or transparency child of a representation whose params
  changed re-runs its `ofScript` query before `areEqual` (`plugin/core/src/state/transforms/structure/effects/transparency.ts:60-79`).
  Skip children whose transform and parent object are unchanged, if the state spec (`core/src/state/_spec`) allows it;
  if it does not, raise it.
- **S8, state interpolation.** `AnimateStateInterpolation` (`state-interpolation.ts:39-71`) compares a cell that holds
  params with `shallowEqual`, which is always false, so every differing transform is updated on every tick; its ref
  match is quadratic.

- **S9, param definitions for dependencies (found by the bench).** `relinkCells` → `getEffectiveDependsOn` →
  `tryGetDefaultParamDefinition` (`core/src/state/transform.ts:189`) calls `params(undefined, undefined)` for each
  scheduled cell, cached on the `Transform`: as with S5, every frame's transforms are new, so the cache never holds. The
  schema without a source depends only on the transformer: cache it there.

## Trajectory playback

- **T1, the model cell is recreated every frame (cause checked; the effect is to be confirmed).** `ModelFromTrajectory`
  (`plugin/core/src/state/transforms/structure/hierarchy.ts:65`) has no `update`, so the default Recreate runs every
  frame. Its `dispose` calls `customProperties.dispose()` on the object that every frame's model shares
  (`model/src/model/structure/model/model.ts:118`, `formats/structure/basic/parser.ts:75`): every step may dispose all
  frames' assets. Confirm it with a test first. Then give it an update that returns Updated with the new frame's model
  and does not dispose shared properties.
- **T2, 05/B1, still open (checked).** `Structure.remapModel` (`model/src/model/structure/structure/structure.ts:421`,
  `parent?.remapModel(m)`) remaps the parent for each component separately. Root-level properties (secondary structure
  and DSSP, ASA through `inheritedPropertyData`) are then computed once per component plus once for the root, every
  frame. The separate roots also have distinct conformations, so `areRootsEquivalent` (`symmetry.ts:170`) may fail across
  components: highlight and selection across components may break after the first frame. Check it with a test. Share
  one remapped root per (parent, model). See [05-model-structure.md](../perf-audit/05-model-structure.md) and
  [bugs.md](../perf-audit/bugs.md), 05/B1.
- **T3, recreated instead of remapped on a coordinate-only change.**
  - `CustomStructureProperties` (`hierarchy.ts:354-366`, `a.data !== b.data`) re-attaches its properties.
  - The overpaint, transparency and other effect cells recreate because `newStructure !== info.structure`
    (`effects/overpaint.ts:67`) and re-run their `ofScript` query for each layer.
  - Expression components re-run their query (`plugin/core/src/state/helpers/structure-component.ts:211-221`).
- **T4, topology caches lost on remap.**
  - `Structure.remapModel` builds new units and a new structure, so `unitSymmetryGroups`, `hashCode`, `models` and
    `carbohydrates` are recomputed although they depend only on topology. Carry them over.
  - `Unit.remapModel` keeps bonds and rings but loses `transientCache` (valence model, interaction features). Check
    which of those read coordinates before carrying anything.
- **T5, behaviors.**
  - `SelectLoci` (`plugin/core/src/behavior/dynamic/representation.ts:209-226`) deselects and re-selects on every
    representation every frame, because its equivalence test fails on new conformations: use topology equivalence.
  - With a focus active, the focus manager (`state/manager/structure/focus.ts:196-215`) makes
    `StructureFocusRepresentation` commit a second state update per frame. A new expression object recreates
    `StructureSelectionFromExpression` (`selection.ts:64`).
  - With the frame index unchanged, an empty update still goes through the task manager on every tick
    (`state/animation/trajectory.ts:137-232`) and toggles `isUpdating`.
  - Prism's marks (`prism/viewer/src/sources/marks.ts:159-201`) rebuild their overpaint once per animation frame. Keep
    it if correctness needs it.
- **T6, visuals (repr and geo layers only; `gl/` is frozen).** Spheres re-upload their group data, and meshes their
  index and group buffers (`graphics/src/geo/geometry/spheres/spheres.ts:98-110`, `mesh/mesh.ts:162-178`), although only
  positions and normals changed. `ComplexVisual` (`graphics/src/repr/structure/complex-visual.ts:182-215`) forces color
  and size updates every frame, because its equivalence test fails on new conformations.

## Both

- **R1, representations rebuilt around unchanged params.** The representation transformer's update
  (`plugin/core/src/state/transforms/structure/representation.ts:129-157`) always releases and re-creates the color and
  size themes, and `resolveParams` re-runs `getApplicableTypes` over every representation, color and size theme
  (:81-99). The MVS multilayer color theme re-runs every layer's selector, not only the changed one
  (`mvs/runtime/src/components/multilayer-color-theme.ts:150-198`). `UnitsRepresentation` re-runs
  `computeChildAwareTransformGroups` without a cache (`graphics/src/repr/structure/units-representation.ts:193`). Skip
  the theme work when its params and data context are unchanged, and cache the params definition by structure topology.
- **Open, not in this plan.**
  - Per-visual `Clip.getClip` allocations, and `ValueCell.update` where `updateIfChanged` would do
    (`graphics/src/geo/geometry/base.ts:247-255`).
  - The scene's bounding sphere recomputed on non-geometric updates (`gl/scene.ts:440-460`: frozen).
  - A transform-matrix frame recreates `TransformStructureConformation`, which has no update
    (`hierarchy.ts:165-183`), so the whole subtree re-runs where instance matrices would do. It is a candidate after
    step 3.
  - MVS step 5, frames computed at playback.

## Questions

- **S6: keep it?** The exosome profile does not show `fromJSON` or the three walks among its top entries; the cost is in
  the evaluated subtrees. Proposal: decide after steps 1–5, by measuring `fromJSON` + `findDeletes` +
  `findUpdateRoots` + `syncNewStates` per frame; drop S6 if they stay under ~10% of a frame. If it is kept, one
  semantic question for `core/src/state/_spec/reconcile.md` §6: `setSnapshot` (a tree passed in) resets every cell's
  state flags to the frame's, while a builder keeps the cells' flags. Applied as a builder, a component the user hides
  during playback stays hidden; today the next frame shows it again. Which one is wanted?
- **S7: in the reconciler, or in the transformers?** `reconcile.md` §5 evaluates every child of an evaluated node.
  Skipping a child whose version equals its cell's when its source returned Unchanged would be safe. But it does not
  cover the case in the finding: the transparency child of a representation updated in place has a changed source, so
  it must run. Proposal: no reconciler change. Instead, the effect transformers' `update` (overpaint, transparency,
  clipping, substance, emissive) skips `ofScript` when the layers' params and the structure are those it last ran on.
  That needs no spec change. Agree?

## How

- **Bench first (step 0).** A headless playback bench in `tests/src/pipeline` (extend the existing one).
  - Measure ms and bytes allocated per frame, as the median and maximum over a played range.
  - Snapshots: the MVS fixtures of the MVS research (motm1, exosome; local only), plus a plain multi-snapshot state
    with a camera transition.
  - Trajectories: an NMR ensemble from `data/` or the test fixtures, and a synthetic multi-frame trajectory of
    ~10k–100k atoms.
  - Representations: cartoon, ball-and-stick and a surface, with and without components, focus and marks.
  - Record the baseline in a table here, and before → after per commit in the PR.
- **One commit per id**, a test that pins today's per-frame output first, then the change. Pinned output means the
  states, params and camera per frame, the number of applies, geometry checksums, and markers. Tests follow the root
  `AGENTS.md` lanes: big fixtures go in medium.
- **The rules:**
  - Run what you touched (`pnpm test:unit <dirs>`, `pnpm test:medium` where it applies), plus `pnpm lint` and
    `npx tsc -b`.
  - Add one entry file in `.v6/changelogs/`.
  - Update the core state and model `_spec`s where semantics move.
  - `gl/` and Canvas3D internals are frozen for the WebGPU work.
- **Review gates:** a fresh reviewer subagent after step 6 and again after step 7, before the PR. It reviews the diff
  against master for per-frame state identity, custom property lifetimes, sharing the root across components, marking
  after many frames, camera transitions, and the fallbacks of the diff-apply.
- **This file is the record:** each step's row moves to done, with its numbers, as it lands. Once the PR is open, the
  rows get the PR number.
