# Slice 10: `prism/viewer/src`

The Prism viewer: the API operations and their run path, the history's checkpoint port, the sources (scene, style
layers, interaction, marks, sequence, state tree), the journal record and the Seq\* integration (`sequence/`). Measured
on `v6-opts-261019` (`90d717f49`), Node, Apple M4, headless (`createPrismViewerHeadless`, no rendering), N × `1cbs`
(1401 cells at N=100) and the 186k-atom `1tup-2x2x2`; method in [00-profile.md](00-profile.md). A second pass re-read
each finding with its callers and re-ran the benches. MVS loading is in `.v6/mvs-optimizations/research/` and is not
repeated. Nothing here touches Canvas3D or `graphics/src/gl`.

## Hot paths

| Path                                                     | When                                                   | Cost (N=100 × 1cbs unless noted)                                                                 |
| -------------------------------------------------------- | ------------------------------------------------------ | ------------------------------------------------------------------------------------------------ |
| `HistoryManager.#take` → `checkpoint()` + `checkpointSize` | before every undoable run (loads included)             | 0.53 + 4.45 ms of a 9.0 ms single-structure `representation.color`; 17 ms at N=300 (P1)          |
| plugin `StructureInfo` on `object.created`               | every created cell (layers, representations, loads)    | first color over all: 16 ms at N=25, 100–111 ms at N=100; 23 % of a load at N≈100 (P4)           |
| plugin measurement `sync`                                | every non-transaction state change                     | ~2 ms per op, 16.5 % of the op profile (P5)                                                      |
| sequence `update` → `build` → link reattach              | every hierarchy emission with a new Structure (frames) | +2.2–2.8 ms per frame step on 11.5k atoms, headless (P3)                                         |
| `readStyleLayers` (memo hit)                             | per state change while the layers list is shown        | 1.12 ms at 300 representations (P6)                                                              |
| `picked` (`pickedOf`)                                    | per state change while Superpose is open               | 3.96 ms for an all-atom selection at 186k atoms (P7)                                             |

Already fine: a no-change scene sync 0.2 ms at 100 structures; `watchedState` reads coalesced per frame; `loadedNames`,
measurement descriptions, style-layer output memoized; lazy contents rows, `contentsLoci` 0.10–0.14 ms and
`contentsSummary` 0.06 ms at 186k atoms; `chainsOf` cached per Structure; `selectionSummary` 0.02 ms; record and History
model windowed; bulk loads in one `dataTransaction`.

## Performance findings

| Id  | Title                                                                         | Where                                                                                                                                   | Kind        | Impact | Estimate                                                                   | Effort | Verdict   | Frozen |
| --- | ----------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ----------- | ------ | -------------------------------------------------------------------------- | ------ | --------- | ------ |
| P1  | Checkpoint size = full `JSON.stringify` with a replacer on every undoable step | [checkpoint.ts:81](../../packages/prism/viewer/src/sources/checkpoint.ts#L81)                                                           | time        | high   | checkpoint + size 5.1 → 1.7 ms (N=100), 17.1 → 4.5 ms (N=300); ~35–40 % of each op | S      | confirmed | no     |
| P4  | (plugin) `StructureInfo` runs 3 `rootsOfType` walks per created object        | [structure-info.ts:133](../../packages/plugin/core/src/behavior/dynamic/custom-props/structure-info.ts#L133)                              | time        | high   | first overpaint/transparency over all at N=100–200: ~3–9× faster; loads ~20 % | S      | confirmed | no     |
| P2  | Each checkpoint re-serializes every transform (~147 KB per step at 1401 cells) | [checkpoint.ts:43](../../packages/prism/viewer/src/sources/checkpoint.ts#L43)                                                           | memory      | medium | 14.7 MB → ~1–2 MB for a full 100-step history at N=100; less GC per op     | S      | confirmed | no     |
| P3  | Every frame step rebuilds the Seq\* document and reattaches the 3D link       | [sequence.ts:152](../../packages/prism/viewer/src/sources/sequence.ts#L152)                                                             | time        | medium | −2.2–2.8 ms/step headless on 11.5k atoms, more with rendering and big structures; no band flash | M      | confirmed | no     |
| P5  | (plugin) measurement `sync` walks the whole tree 6× per state change          | [measurement.ts:435](../../packages/plugin/core/src/state/manager/structure/measurement.ts#L435)                                         | time        | medium | ~2 ms (≈17–25 %) off every state-changing op at 1401 cells, O(cells)       | S      | confirmed | no     |
| P6  | `readStyleLayers` runs 5 selects per representation even on a memo hit        | [style-layers.ts:145](../../packages/prism/viewer/src/sources/style-layers.ts#L145)                                                     | time        | low    | 1.1 → ~0.1 ms per change at 300 representations                            | S      | confirmed | no     |
| P7  | Superpose `picked` recomputes `pickedOf` on every state change                | [interaction.ts:64](../../packages/prism/viewer/src/sources/interaction.ts#L64)                                                         | time        | low    | ~4 ms per frame on a 186k-atom selection while the section is open         | S      | confirmed | no     |
| P8  | `trajectory.frame` seeks are written to the record, one per seek              | [trajectory.ts:61](../../packages/prism/viewer/src/api/trajectory.ts#L61)                                                               | time+memory | low    | hundreds fewer record writes (IndexedDB when kept) per scrub               | S      | confirmed | no     |
| P9  | `chainsOf` indexes heavy atoms of water chains nothing links                  | [chains.ts:93](../../packages/prism/viewer/src/sequence/molstar/chains.ts#L93)                                                          | time+memory | low    | 1.3–1.4 → 0.85 ms per new Structure at 186k atoms (12k waters)             | S      | confirmed | no     |
| P11 | Link's 3D `onSelect` reacts per structure event, not once a task              | [sequence.ts:285](../../packages/prism/viewer/src/sources/sequence.ts#L285)                                                             | time        | low    | S remaps of the target selection per bulk select → 1                      | S      | confirmed | no     |
| P10 | `StateTreeSource.node()/details()` rebuild the tree on each call when unwatched | [source.ts:318](../../packages/prism/viewer/src/features/state-tree/source.ts#L318)                                                     | time        | low    | k unwatched reads: k × O(cells) → 1 per tree version; no hot caller found  | S      | plausible | no     |

### P1: checkpoint size by full JSON.stringify

- **Mechanism**: `#take` ([history.ts:608](../../packages/prism/ui/src/managers/history.ts#L608)) calls
  `port.checkpoint()`, times it, then `port.size()` = `checkpointSize`: `JSON.stringify` of the whole
  `PluginState.Snapshot` with a replacer called per value. O(cells) on every undoable run (each load too), so N
  one-at-a-time loads are O(N²) in total. The time budget (`ms`, 50/500) covers only `checkpoint()`, so the size pass is
  invisible to the "slow checkpoints turn undo off" rule. Params are immutable and shared by reference between snapshots
  (`StateTransform.toJSON` uses `_id`), yet re-serialized each step.
- **Evidence**: `bench-op100.ts`, N=100: checkpoint 0.53 ms, size 4.45 ms (1.07 MB); `representation.color` on one
  structure 9.0 ms; profile `checkpointSize` 30.8 % of the run (≈4.5 ms per op, half the op). `bench-cksize.ts`: per-params
  memo prototype 5.09 → 1.68 ms at N=100 (1,066,607 vs 1,051,165 B, −1.4 %), 17.1 → 4.5 ms at N=300. Load profile at
  N≈100: 8.5 %. The 5 MB budget turns undo off at ~470 structures of this size.
- **Fix**: `WeakMap<params, number>` of `JSON.stringify(params, replacer).length`; sum it with the small fixed fields
  (ref, parent, transformer, version, tags, state) and stringify only `behavior`/`canvas3d`/`custom`/non-tree parts each
  time. Better: a running total kept from cell created/removed/updated events. Count the size pass in `#take`'s `ms`, or
  say in sequence.md §7.4 that it is excluded.
- **Risk**: budget bytes move by ~1.5 %; byte-size assertions in prism-ui history and checkpoint tests may need a new
  baseline. `customSerialization` params are fresh objects each snapshot: key those on the transform instead.

### P4: StructureInfo selects on every created object (plugin, cross-slice)

- **Mechanism**: on every `events.object.created` of any type it runs `handle()`, then `setModelMaxIndex()` (inside it
  `maxModelIndex`: a second Model walk) and `setStructureMaxIndex()`. Correction to the audit: `rootsOfType` stops
  descending at the first match, so a walk visits the nodes above models/structures (~5 per loaded structure), not every
  cell, through immutable.js child iteration. Still O(N) per created object, O(K·N) for K creations; an "overpaint all"
  creates one layer cell per representation (K ≈ 3N), so O(N²).
- **Evidence**: `bench-colorall.ts` re-run: first `representation.color` with a color over all: 16 ms at N=25,
  100–111 ms at N=100 (auditor: 372 ms at N=200); `representation.transparency` 10 / 78–83 ms. Superlinear: 4× N gives
  6–8× time. Profile at N=100: structure-info.ts:133 27.4 % of the run (≈55 % of busy time), `setStructureMaxIndex`
  15.5 %, `maxModelIndex` 11.8 %. Load profile at N≈100: 23.4 %.
- **Fix**: call `setModelMaxIndex`/`setStructureMaxIndex` only when the created object is a non-decorator `Model` /
  root `Structure` (the only ones that can move the max); keep the max indices in fields, updated on create and remove.
  `handleModel`/`handleStructure` already assign the index.
- **Risk**: plugin/core owner. Keep `MaxIndex` semantics on removals (today a removal does not lower it until the next
  create; update on `object.removed` if a test relies on lowering).

### P2: checkpoints duplicate unchanged transforms

- **Mechanism**: `checkpoint()` → `plugin.state.getSnapshot` → `StateTree.toJSON` allocates a new `Serialized` per
  transform (plus its `state` object) every step; only params are shared. Capacity 100 plus redo checkpoints.
- **Evidence**: `bench-ckmem.ts --expose-gc`: 100 checkpoints retain 14.7 MB at N=100 (147 KB each); 3.6 MB at N=25.
  Small next to a 100-structure scene's own heap, but it is GC churn on every op.
- **Fix**: in the port, `Map<ref, {version, serialized}>` from the previous checkpoint; swap in the cached object when
  `(ref, version)` is unchanged (a new version is a new transform). Or memoize `StateTransform.toJSON` per transform
  (`WeakMap`) in core's `StateTree.toJSON` (core-state slice, not frozen). Share the key with P1's size cache.
- **Risk**: consumers must not mutate serialized transforms: `keepFrames` copies before changing; `setSnapshot` →
  `StateTransform.fromJSON` reads only. Undo/redo tests cover it.

### P3: sequence rebuild on every frame step

- **Mechanism**: `update()` keys `HostStructure` by `keyOf(ref, data)`, `data` the Structure object; a frame step makes
  a new Structure, so every hierarchy emission gives a new key: a new `StructureSource` (empty `structures` cache),
  `Viewer.build()` ([model.ts:659](../../packages/prism/viewer/src/sequence/model.ts#L659)): `sequenceDocument` +
  `seq.setDocument` (complex provider, `chainsOf`, `complexInputs`, tracks, behaviors' patches), `structureState` set to
  undefined then the new value (:694–703, the band empties and refills), and `linkSequence`
  ([link.ts:57](../../packages/prism/viewer/src/sequence/link.ts#L57)) drops and reattaches `StructureLink`. The plugin's
  model-index animation goes the same way.
- **Evidence**: `bench-frames.ts`, 10-model 7qpd (11.5k atoms/model), 20 `trajectory.frame by 1`: 20 source emissions,
  40 viewer structure changes; median step 14.0 ms without sequence, 16.2 ms with, 16.8 ms with the link. (The audit's
  168 system registrations did not reproduce with this hook; the reattach is evident from the code.) Seq\* relayout and
  repaint come on top in the page; `chainsOf` alone is 1.3–2.4 ms per new 186k-atom Structure.
- **Fix**: key by ref + trajectory + a hierarchy-equality check (`Model.areHierarchiesEqual` with the previous model,
  or `Structure.areUnitIdsAndIndicesEqual`); pass the current Structure separately. In `linkSequence`, when only the
  Structure object changed, rebind `StructureLink` to it (the view3d `onto` remap compares `root`, so the link must
  follow the new object) instead of drop + attach. One `StructureSource` per Trajectory (`WeakMap`, as `data-by-url`
  does).
- **Risk**: sequence.md §2 "A new key rebuilds the view" needs a same-hierarchy rule. Tests: frame steps keep the
  document; a model with another hierarchy rebuilds; hover/selection map onto the frame shown.

### P5: measurement sync (plugin, cross-slice)

- **Mechanism**: `sync()` on every non-transaction `events.changed` calls `getTransforms` 6×; each starts with
  `findTagInSubtree(tree, Root, MeasurementGroupTag)`, which visits every node when no group exists (the usual case).
- **Evidence**: `bench-op100.ts` profile: `sync` 16.5 % of the run, `findTagInSubtree` 15.5 % (≈2 ms of a 9 ms op,
  profiler-inflated; the audit's 1.5 ms is the floor).
- **Fix**: find `groupRef` once per sync, cached by `state.tree` identity (immutable); no group → the six empty lists.
- **Risk**: plugin/core owner; none to semantics.

### Low

- **P6**: per representation × 5 kinds, `StateSelection.Generators.ofTransformer(t, repr).withTag(tag)` is compiled and
  run before the memo compares transforms. `bench-layers.ts`: 1.12 ms memo hit, 1.62 ms cold at 300 representations.
  Fix: one pass over `tree.children.get(repr)` matching the 5 transformers through a Map; or invalidate the found list
  only on cell events of those transformers. Keep the order (kind, then as found).
- **P7**: `picked` follows `merge(selection.events.changed, state.data.events.changed)`; `pickedOf` is O(selected atoms)
  and allocates. `bench-struct.ts`: 3.96 ms on an all-atom selection at 186k atoms. The `selection` source beside it
  already avoids data events for this reason. Fix: follow `selection.events.changed` + `hierarchy.behaviors.selection`,
  or memoize per (structure ref, loci object).
- **P8**: `trajectory.frame` is `history: false` but not `record: false`; `isQuiet`
  ([record.ts:216](../../packages/prism/viewer/src/journal/record.ts#L216)) misses `trajectory.*`. The scrubber
  ([ui/trajectory.ts:51](../../packages/prism/viewer/src/ui/trajectory.ts#L51)) and scene seek
  ([scene.ts:1273](../../packages/prism/viewer/src/sources/scene.ts#L1273)) run it once per completed seek, so a drag
  writes ~one entry per frame (`summarizeValues`, append, listeners; IndexedDB when kept). journal.md lists `record:
  false` as "pointer-rate view state". Fix: coalesce consecutive same-origin `trajectory.frame` runs into one entry (the
  agent keeps the frame shown), or `record: false` (as `snapshot.seek`). Policy call.
- **P9**: non-polymer branch allocates a `Set` per residue and compares names for waters too; both consumers skip
  `kind === 'water'` ([link.ts:117](../../packages/prism/viewer/src/sequence/molstar/link.ts#L117),
  [complex.ts:201](../../packages/prism/viewer/src/sequence/molstar/complex.ts#L201)); `seq-units` looks up only chains
  that became units. `bench-chains-water.ts`: 1.27–1.44 → 0.85–0.88 ms on 1tup-2x2x2. Fix: water chains with count 0,
  `chainAt`/`elementAt` left −1.
- **P11**: `onSelect` subscribes to `selection.events.changed` directly; `api/select.ts` `write` emits once per edited
  structure, so a bulk select over S structures remaps the target selection and runs `Loci.areEqual` S times. The
  source's own listener (:172) uses `onceATask`. Fix: wrap in `onceATask`.
- **P10**: `current()` sets `all = true` and flushes whenever nothing listens, so each unwatched `node()`/`details()`/
  `roots.value` is an O(cells) rebuild. Callers found: `state.get` (once per agent call) and the State tab model, whose
  calls happen while the tab watches. Fix: remember the `state.tree` (and a cell-object version) the cache was built
  from; rebuild only when it changed.

## Bugs

| Id  | Title                                                                                 | Where                                                                                       | Severity | Verdict   |
| --- | ------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | -------- | --------- |
| B1  | Followed structure's source has `frame: null`: no SIFTS numberings, no 'predicted' origin, no scheme auth labels | [sequence.ts:163](../../packages/prism/viewer/src/sources/sequence.ts#L163)                 | medium   | confirmed |
| B2  | `paintMarks` keeps removed Representations until the next paint                       | [marks.ts:143](../../packages/prism/viewer/src/sources/marks.ts#L143)                       | low      | confirmed |
| B3  | Whole-chain fallback spans untouched chains between touched ones in a merged unit     | [link-views.ts:94](../../packages/prism/viewer/src/sequence/molstar/link-views.ts#L94)      | low      | confirmed |

- **B1**: `update()` builds `{ format, trajectory, frame: null, structures: new Map() }`; the loaders' `parsed` hook
  ([model.ts:383](../../packages/prism/viewer/src/sequence/model.ts#L383)) hands it to the `molstar:structure` complex,
  and `complexInputs` reads `frame` for `ma_qa_metric_local` (origin), `siftsNumberings` (UniProt numberings, default
  'best') and `pdbx_poly_seq_scheme` (auth labels of unobserved residues)
  ([complex.ts:196](../../packages/prism/viewer/src/sequence/molstar/complex.ts#L196)). All silently dropped.
  `data-by-url.ts:110–113` passes the CIF block, so only the main band is affected. Failure (`check-sifts2.ts`):
  `1tup.bcif` molecule 3 has 1 UniProt numbering with the frame, 0 with null; `AF-P04637-F1-model_v6.bcif` origin
  'predicted' with the frame, 'experimental' with null; `check-sifts.ts` headless: the assembled set has no
  `uniprot:P04637` string. Effect: no UniProt numbering/mapping/alignment and no predicted-model styling on the
  followed structure. Fix: frame from the CIF block above the trajectory cell (`SO.Format.Cif` parent → `blocks[0]`),
  or reuse `dataByUrl.structure(address)`; cache the source per Trajectory (also serves P3). Test: uniprot numbering on
  1tup, 'predicted' on the AF fixture. Some sequence snapshots will change (intended).
- **B2**: `painted`/`unpainted` are replaced only in `paint()`, scheduled on marks/colors changes and on
  `object.created/updated` of Representation3D(State). Removing a representation of a marked structure fires none, so
  the removed `Representation` stays referenced until the next paint (or the marks clear). Correction to the audit:
  `destroy()` clears its visuals ([units-representation.ts:397](../../packages/graphics/src/repr/structure/units-representation.ts#L397)),
  so meshes and render objects are not kept; the closure keeps its Structure, theme and props. Bounded, temporary.
  Fix: subscribe to `cell.removed` and drop entries by ref, or store refs and resolve the repr at paint time.
- **B3**: `wholeChains` keeps one `[lo, hi)` per Mol\* unit over all touched chains. Mol\* merges consecutive one-atom
  chains (ions, linked as 'ligand') into one unit (chains.ts:10). Failure: a Seq\* hover over 2000 items (or a select
  over 20000) that touches ions A and C of one unit but not B lights/selects B too. Fix: one interval per touched chain,
  the unit's indices as their sorted union.

## Checked and dropped

- Nothing refuted. Corrections: P4's walks stop at the matched type (O(structures), not O(cells), per created object;
  still O(N²) for "all" ops); P5 measured nearer 2 ms than 1.5 ms; P3's 168 registrations not reproduced (frame cost and
  40 structure changes are); B2 does not retain meshes (destroy clears visuals), downgraded to low; P9 measured
  (−35 %); P10 left plausible (no hot caller).

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/prism-viewer/`,
bundled with `node ../build.mjs X.ts X.mjs` (esbuild, `molstar-src` condition), run with `node X.mjs [N]`; profiles
read with `node prof.mjs X.cpuprofile X.mjs.map`.

| Script                  | What                                                                                     |
| ----------------------- | ---------------------------------------------------------------------------------------- |
| `bench-op100.ts`        | P1, P5: checkpoint/size split and a profiled single-structure `representation.color`     |
| `bench-cksize.ts`       | P1: `checkpointSize` vs per-params memo prototype (N=100, 300)                           |
| `bench-ckmem.ts`        | P2: heap retained by 100 checkpoints (`node --expose-gc`)                                |
| `bench-colorall.ts`     | P4: first overpaint/transparency over all structures, `prof` for a profile               |
| `bench-load100.ts`      | P4, P1: profiled one-at-a-time loads at N≈100                                            |
| `bench-frames.ts`       | P3: 20 frame steps on `multi.cif` (10-model 7qpd), `noseq` / default / `seq link`        |
| `check-sifts2.ts`       | B1: `complexInputs` with and without the CIF frame on 1tup and AF-P04637                 |
| `check-sifts.ts`        | B1: headless viewer, assembled set has no `uniprot:P04637`                               |
| `bench-layers.ts`       | P6: `readStyleLayers` memo hit and cold at 300 representations                           |
| `bench-struct.ts`       | P7, P9, hot-path baselines: `chainsOf`, `pickedOf`, `summarize`, contents rows at 186k    |
| `bench-chains-water.ts` | P9: `chainsOf` copy (`chains-v.ts`) with and without water indexing                      |
| `bench-scene.ts`, `bench-ops.ts`, `bench-seqbuild.ts` | baselines: scene sync, op costs, Seq\* document build             |
