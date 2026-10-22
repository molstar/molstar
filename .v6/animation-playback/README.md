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
| 1    | S1–S5, S9              | Snapshot playback: the per-tick re-apply bug and per-tick waste           | open   |
| 2    | T1, T2                 | Trajectory: the model cell and the shared root (correctness first)        | open   |
| 3    | T3, T4                 | Trajectory: remap instead of recreate; keep the topology caches           | open   |
| 4    | R1                     | Representations: no theme or params rebuild when nothing they read moved  | open   |
| 5    | T5                     | Behaviors: selection, focus, the empty update, Prism marks                | open   |
| 6    | T6                     | Visuals: unchanged buffers stay, complex visuals keep colors              | open   |
| 7    | S6, S7                 | A frame applied as a diff of the changed refs; unchanged children skipped | open   |
| 8    | S8                     | State interpolation: equality and the ref match                           | open   |

Order: the bench, then the bugs (1, 2), then the work per frame from the bottom of the tree to the top (3–6), then the
state-level change (7). Step 7 changes how a snapshot is applied, so it needs the core state spec, and the steps before
it make its gain measurable on its own.

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
