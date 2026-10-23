# Slice 10: `prism/viewer/src` (v2)

The Prism viewer: API operations and their run path, the history's checkpoint port, the sources (scene, style layers,
interaction, marks, sequence, state tree, clipping, structure move), the journal record and the Seq\* integration
(`sequence/`). Measured at `0944823d4`, Node 25.6, Apple M4, headless (`createPrismViewerHeadless`, no rendering),
N × `1cbs` and the 10-model `7qpd` (`multi.cif`). MVS loading/animation is in `.v6/mvs-optimizations/research/` and is
not repeated. P12 sits in frozen Canvas3D/`graphics/src/gl`.

Since v1: #2008 fixed nearly all of v1 (checkpoints, frame-step rebind, style layers, Superpose memo,
`trajectory.frame` record, water chains, link `onSelect`, `frame: null`, marks prune, whole chains); #2013 finished
StructureInfo (P4) and the measurement group cache (P5, partly); #2009 added the Log's Transforms list; #2021 added
clipping, gizmos and structure move; #2022 animations with `doNotLogTiming`.

## Hot paths now

| Path | When | Cost (headless) |
| --- | --- | --- |
| `Canvas3D.update(repr)` → full precompile record (frozen) | per representation: move-gizmo drag (each pointer move), marks paint, plugin `SyncStructureRepresentation3DState` | ~6.5 µs per render object in the scene per call (11 µs in the bench's schema-texture proxy); drag at N=100 ≈ 12 ms per move; color over all at N=100 ≈ 0.3 s, O(N²) (P12) |
| checkpoint + size | every undoable run | 2.6 + 0.7 ms at 4201 cells; size 0.27 ms at 1401 (was 5.0) |
| single-structure op (`representation.color`, `select.set`) | person/agent | 2.9 ms at N=100, 4.3 ms at N=300; select ~2 ms flat |
| op over all structures | person/agent | linear: color 6 → 21 ms, select 5 → 14 ms for N=50 → 200 |
| load | per structure | 25 ms at N=300; Prism 2-3 % of the profile, rest the plugin's hierarchy walk (14 open item 2) |
| frame step with the sequence linked | playback, scrub | Prism ~2 % of a step; 0 system registrations; but 1 `highlight.clear` (+1 `select.set` with a Seq\* selection) and ~9 log lines per step (P13, B6, P14) |
| plugin measurement `sync` | every non-transaction state change | one `findTagInSubtree` miss per change: 1.15 ms at 4201 cells, 6.9 % of a `focus.set` at N=300 (P5) |
| `findCuts` | per state change while Clipping watches | O(K²) in representations of one cut: 8 ms at K=1500 (P15) |

Cheap at rest: clipping, gizmos and Kinemage sources follow the plugin only while watched; clip previews write uniforms
only.

## v1 → now

Counts: 15 v1 items (P1-P11, B1-B3, 14/NF2 Prism part): fixed 12, partly fixed 1, still open 1, regressed 0, dropped 1.

| v1 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P1 | Checkpoint size = full `JSON.stringify` per undoable step | fixed | #2008 | per-`Serialized` size memo ([checkpoint.ts:136](../../packages/prism/viewer/src/sources/checkpoint.ts#L136)); 5.00 → 0.27 ms at 1401 cells, 13.65 → 0.71 ms at 4201; counted in the history budget |
| P2 | Each checkpoint re-serializes every transform | fixed | #2008 | `WeakMap<StateTransform, Serialized>` shared across checkpoints, flags re-checked per read; 100 checkpoints 14.8 → 1.6 MB; undo capped at 10 |
| P3 | Frame step rebuilds the Seq\* document and reattaches the link | fixed | #2008 | `sameAtoms` keeps the key ([sequence.ts:179](../../packages/prism/viewer/src/sources/sequence.ts#L179)); `StructureLink.rebind`; 20 steps, 0 registrations. Views still restart per step: P13/B6 |
| P4 | (plugin) StructureInfo 3 `rootsOfType` walks per created object | fixed | #2008, #2013, #2026 | roots tracked from events; `rootObjects` 1.0 %, `setMaxIndices` 0.4 % of a load at N=300; #2026: `setMaxIndices` for root models only (12: 8.6 % of many-500 `repr:auto` → 0) |
| P5 | (plugin) measurement sync walks the tree 6× per change | **fixed** | #2008, #2013, #2026 | 6 walks → 1 (#2008); #2026: group ref from cell events, no walk per change (many-500 theme:one 5.2 → 0 ms profiled) |
| P6 | `readStyleLayers` 5 selects per representation on a memo hit | fixed | #2008 | one `doPreOrder` per representation, kinds by transformer Map ([style-layers.ts:127](../../packages/prism/viewer/src/sources/style-layers.ts#L127)) |
| P7 | Superpose `picked` recomputes `pickedOf` per state change | fixed | #2008 | `pickedMemo` per selection loci ([interaction.ts:67](../../packages/prism/viewer/src/sources/interaction.ts#L67)) |
| P8 | `trajectory.frame` seeks written to the record | fixed | #2008 | `record: false` ([trajectory.ts:80](../../packages/prism/viewer/src/api/trajectory.ts#L80)); the link's replay now writes `select.set` per step instead (B6) |
| P9 | `chainsOf` indexes water heavy atoms | fixed | #2008 | water chains with count 0 ([chains.ts:91](../../packages/prism/viewer/src/sequence/molstar/chains.ts#L91)) |
| P10 | `StateTreeSource.node()/details()` rebuild when unwatched | still open | | [source.ts:318](../../packages/prism/viewer/src/features/state-tree/source.ts#L318) unchanged; no hot caller |
| P11 | Link's 3D `onSelect` per structure event | fixed | #2008 | `onceATask` ([sequence.ts:317](../../packages/prism/viewer/src/sources/sequence.ts#L317)) |
| B1 | Followed structure's source has `frame: null` | fixed | #2008 | `cifFrameOf` of the trajectory cell, one source per Trajectory ([sequence.ts:122](../../packages/prism/viewer/src/sources/sequence.ts#L122)) |
| B2 | `paintMarks` keeps removed Representations | fixed | #2008 | `cell.removed` → `prune()` once a microtask ([marks.ts:242](../../packages/prism/viewer/src/sources/marks.ts#L242)) |
| B3 | Whole-chain fallback spans untouched chains in a merged unit | fixed | #2008 | one interval per touched `ChainInstance` ([link-views.ts:99](../../packages/prism/viewer/src/sequence/molstar/link-views.ts#L99)); B5 is a different gap there |
| 14/NF2 | Prism's own state-change subscribers not measured (14) | dropped | | measured: Prism 2-3 % of a load at N=300, per-op costs flat, nothing quadratic on Prism's side |

## Open findings (ranked)

| Rank | Id | vs v1 | Title | Where | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P12 | new | `Canvas3D.update(repr)` re-records every renderable's precompile signature; callers loop per representation | [canvas3d.ts:2240](../../packages/graphics/src/canvas3d/canvas3d.ts#L2240) | time | high | drag N=100 ~12 ms → <0.5 ms per move; color over all O(N²) → O(N): ~0.3 s at N=100, ~8 s at N=500 | S | confirmed | yes |
| 2 | P13 | new | Frame step restarts the link views: `highlight.clear` per frame, `select.set` replay with a Seq\* selection | [link.ts:104](../../packages/prism/viewer/src/sequence/link.ts#L104) | time | low | 1-2 API runs and mark passes per frame fewer | S | confirmed | no |
| 3 | P5 | **fixed** (#2026) | (plugin) measurement sync: one full walk per state change with no group | [measurement.ts:426](../../packages/plugin/core/src/state/manager/structure/measurement.ts#L426) | time | low | ~1.15 ms per change at 4201 cells | S | confirmed | no |
| 4 | P15 | new | `findCuts` gathers targets with a linear `.some`: O(K²) per read | [cuts.ts:101](../../packages/prism/viewer/src/sources/clipping/cuts.ts#L101) | time | low | 8 ms → 0.1 ms at K=1500 per state change while shown | S | confirmed | no |
| 5 | P14 | new; **fixed** #2032 | `trajectory.frame`/`particles.frame` log ~9 transform lines per seek | [trajectory.ts:102](../../packages/prism/viewer/src/api/trajectory.ts#L102) | time+memory | low | 9 log entries per seek fewer; plugin log stops growing with scrubbing | S | confirmed | no |
| 6 | P16 | new | Bulk runs' `did` is one line per structure, kept uncut in the record and status | [record.ts:520](../../packages/prism/viewer/src/journal/record.ts#L520) | memory | low | ~10 KB at N=200 → O(1) per entry | S | confirmed | no |
| 7 | P10 | still open | `StateTreeSource` rebuilds on each unwatched read | [source.ts:318](../../packages/prism/viewer/src/features/state-tree/source.ts#L318) | time | low | k × O(cells) → one per tree version | S | plausible | no |

### P12: Canvas3D.update(repr) records the whole scene (frozen)

- **Mechanism**: `update(repr)` runs `scene.update(repr.renderObjects)`, then `precompileRecorder.invalidate()` (a new
  `WeakSet`) and `record(scene.renderables)`: `getObjectSignature` (every `define` key of the schema, linear textures,
  draw flags) and `getSignatureKey` (sort + joins) for every renderable in the scene, though one representation
  changed. Callers loop over representations, so O(R_changed × R_scene):
  - structure-move preview, per pointer move: `#draw` → `update` per drawn representation
    ([gizmos.ts:173](../../packages/prism/viewer/src/sources/structure-move/gizmos.ts#L173), #2021);
  - `paintMarks`, per painted representation ([marks.ts:201](../../packages/prism/viewer/src/sources/marks.ts#L201));
  - plugin `SyncStructureRepresentation3DState`, per Representation3DState created/updated/removed
    ([representation.ts:53](../../packages/plugin/core/src/behavior/static/representation.ts#L53)).
- **Evidence**: `bench-sig.ts` (signature + key over the scene's render objects): 12 µs/object, 3.05 ms per pass at
  N=50 (250 objects), 5.6 ms at N=100 (500). `representation.color` over all at N=50 updates 50 Representation3DStates →
  50 passes. Correction (`bench-sig2.ts`): the bench's `getRenderObjectSignature` derives linear textures from the schema
  (5.1 µs/object); the canvas's `getObjectSignature` reads `renderItem.getLinearTextures()` instead, so the real pass is
  ~6-7 µs/object: ~3 ms per pass at N=100, a drag move over 4 representations ~12 ms, color over all at N=100 ~0.3 s,
  N=500 (~2500 objects) ~8 s. Not measured in a browser (no native gl here); 14's stub canvas had no recorder.
- **Fix**: frozen. Lazy recorder: `update(repr)` drops only that repr's renderables from `recorded`, or sets a dirty flag,
  and `record()` runs in `getPrecompileList()` (its only reader) or at the next commit. Not frozen interim: the drag
  sets state on all drawn representations and calls `canvas3d.update()` once (one pass; it updates every renderable, so
  measure); the plugin behavior could batch its updates once a task. Test: signatures computed per `update`.
- **Risk**: every drawn signature must still be in the list before `getPrecompileList`; existing precompile tests.

### P13: frame step restarts the sequence link's views

- **Mechanism**: the frames subscription ([link.ts:104](../../packages/prism/viewer/src/sequence/link.ts#L104))
  rebinds `StructureLink`, then `stop()` + `start()`. `stop()` calls `lit?.highlight(null)` unconditionally
  ([link.ts:55](../../packages/prism/viewer/src/sequence/link.ts#L55)) → `link.highlight.clear()` (clears every plugin
  highlight, marks hover). `start()` runs `linkViews` again: its replay (`first=true`, `seen=null`,
  [link-views.ts:183](../../packages/prism/viewer/src/sequence/molstar/link-views.ts#L183)) resends a non-`host:3d`
  selection as `link.select.set`. Also a new `SequenceView3D` (3 channels) and a label provider add/remove per frame.
- **Evidence**: `bench-frames2.ts` (7qpd, 20 `trajectory.frame by 1`): `highlight.clear/link` 20; with a Seq\*-origin
  selection also `select.set/link` 20 and 20 record entries; 0 registrations; selection-changed 0 (the runs change
  nothing). Headless Prism ~2 % of a step; in the page each run adds a Canvas3D mark pass, small for one trajectory
  structure. Impact lowered from the auditor's medium to low: the visible effects are B6.
- **Fix**: on a rebind keep the views and swap the structure `onto`/`fromLoci` use (a `rebind(structure)` on the views),
  or `start(…, { rebind: true })`: no replay, and `stop()` clears only a highlight this link set (track non-null).
- **Risk**: keep the replay for a real attach and the clear when the followed structure changes.

### Low

- **P5** (plugin): `groupRef()` caches by `tree.transforms`, new on every edit, so a miss walks every node
  (`findTagInSubtree` never stops early, 03/B10): 1.15 ms at 4201 cells, 6.9 % of `focus.set` at N=300. The comment at
  [measurement.ts:448](../../packages/plugin/core/src/state/manager/structure/measurement.ts#L448) ("no walk of the
  tree at all") is wrong. Fix: keep the cached ref while `tree.transforms.get(ref)` still has the tag; with none, rescan
  only when a cell tagged `MeasurementGroupTag` was created (an `events.cell.created` flag). Correct the comment.
  **Fixed (#2026):** the ref from tagged `cell.created`/`cell.stateUpdated`, checked O(1) per sync, one search after
  the group goes or loses its tag; comment rewritten. many-500 theme:one 58 → 46 ms, remove:one 52 → 37 (medians, with
  12/10/P4).
- **P15**: per clippable cell, `cuts.find` (O(cuts)) and `cut.targets.some` (O(targets)); `clipping.add` with
  `everything` puts one object on every representation (`allCells`), so K = all representations. `bench-cuts.mjs`:
  0.27 / 8.1 / 15.3 ms at K=300 / 1500 / 3000 vs 0.04-0.3 ms with a `Set`. `ClippingApi.cuts` re-reads per state change
  while watched. Fix: a `Set` of refs per cut, a Map by cut key; memoize `keptBy`. Keep tree order (clipping.md §1).
- **P14**: seeks call `plugin.state.updateTransform` without `doNotLogTiming`; core logs "Updated X in N ms" per cell,
  `plugin.log.entries` is never trimmed ([context.ts:692](../../packages/plugin/core/src/context.ts#L692)) and
  [viewer.ts:499](../../packages/prism/viewer/src/viewer.ts#L499) parses each into the Log's Transforms (keep 100).
  `bench-frames2`: 180 log events for 20 steps; ~11 seeks push out every real transform line. Fix:
  `PluginCommands.State.Update(…, { options: { doNotLogTiming: true } })` in `trajectory.frame`, `particles.frame`
  ([operations.ts:384](../../packages/prism/viewer/src/features/particles/operations.ts#L384)), one helper next to
  `whileViewing`.
- **P16**: `describe` ([select.ts:351](../../packages/prism/viewer/src/api/select.ts#L351)) joins one summary per
  structure; the record keeps `did` uncut (values are summarized) and kept sessions write it to IndexedDB; the status
  posts it whole. ~10 KB at N=200 (`bench-ops`). Fix: condense past 3 structures, `plainText(did, ~1000)` in the
  record; keep per-structure detail in `value` for the agent.
- **P10**: `current()` sets `all = true` and flushes whenever nothing listens: each unwatched `node()`/`details()`/
  `roots.value` is O(cells). Callers: `state.get` once per agent call; the State tab reads while watching. Fix: remember
  the `state.tree` and an object-version counter; rebuild on change. Plausible: no hot caller.

## Bugs

| Id | vs v1 | Severity | Title | Where | Verdict |
| --- | --- | --- | --- | --- | --- |
| B4 | new | medium | While any animation (camera spin included) or a snapshot plays, outside changes never invalidate checkpoints | [checkpoint.ts:205](../../packages/prism/viewer/src/sources/checkpoint.ts#L205) | confirmed (code) |
| B6 | new | low | Linked sequence: each frame step clears the 3D highlight and, with a Seq\* selection, writes a `select.set` record entry | [link-views.ts:183](../../packages/prism/viewer/src/sequence/molstar/link-views.ts#L183) | confirmed |
| B5 | new; **fixed** #2034 | low | `wholeChains` drops a non-polymer chain's hydrogens after its last heavy atom | [link-views.ts:107](../../packages/prism/viewer/src/sequence/molstar/link-views.ts#L107) | confirmed (code) |

- **B4**: `watchOutsideChanges` returns on `isViewing(plugin)` = the viewing counter or
  `managers.animation.isAnimating` ([trajectory.ts:176](../../packages/prism/viewer/src/sources/trajectory.ts#L176)),
  and on `snapshot.state.isPlaying`: exemption by time, not origin. `animation.play` offers `animate-camera-spin` and
  rock, which change no data and run indefinitely. Failure: Camera Spin plays; an extension or page script commits
  `plugin.build().to(…).apply(…)`; no invalidate; Undo restores the pre-spin checkpoint and the outside edit vanishes
  without the "can no longer be undone" notice (sequence.md §7.3). Same, briefly, for any outside commit during a
  `whileViewing` seek's await. Fix: exempt the animation's own updates (frame-changing animations by name, or tag the
  refs they write); camera animations exempt nothing. Test: outside commit during camera spin invalidates; model-index
  playback does not.
- **B6**: P13's restart. `highlight.clear` → `lociHighlights.clearHighlights()` drops the person's 3D hover highlight each
  step (it may also drop on the visual update itself; not checked in a page). The replay sends a `pointer`/`keyboard`
  Seq\* selection as `select.set` (origin link, recorded: not `record: false`), so 20 steps → 20 record entries
  (`FROMSEQ=1 bench-frames2`), against P8's intent and journal.md's pointer-rate rule. Fix as P13; test: Seq\* selection +
  20 steps → 0 entries, a 3D hover survives.
- **B5**: `chainsOf` indexes non-polymer chains by heavy atom; `wholeChains` takes `[start[0], end[count-1])`, from the
  first to just past the last heavy atom ([chains.ts:94](../../packages/prism/viewer/src/sequence/molstar/chains.ts#L94)).
  wwPDB lists a ligand residue's hydrogens after its heavy atoms. Failure: a Seq\* select over 20000 items (hover over
  2000) touching a hydrogen-bearing ligand or glycan chain leaves the last residue's trailing hydrogens unselected, while
  the within-limits path (`wholeResidue`) includes them. Fix: span from the chain's atom segment in the unit (as polymer
  residues map), or extend the result to whole residues. Water stays excluded.

## Checked and dropped

- Nothing refuted. Corrections: P12's per-object cost ~6-7 µs, not 11 (the bench's schema-derived linear textures cost
  5 µs/object the canvas does not pay); estimates halved, still O(N²). P13 lowered to low (its visible effects are B6).
  B6's hover loss confirmed from code only.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/prism-viewer/`,
bundled with `node ../build-pv.mjs X.ts X.mjs` (esbuild, `molstar-src`), run `node X.mjs [N]`; profiles with
`node prof.mjs X.cpuprofile X.mjs.map`.

| Script | What |
| --- | --- |
| `bench-ck.ts` | P1, P2, P5: checkpoint/size vs full stringify at N=100/300; `findTagInSubtree` miss |
| `bench-ops.ts` | per-op costs vs N, load profile, P16 `did` sizes |
| `bench-frames2.ts` | P3, P13, B6, P14: 20 frame steps on `multi.cif`; `NOLINK=1`, `FROMSEQ=1` |
| `bench-frameprof.ts` | frame-step profile (Prism share) |
| `bench-sig.ts` | P12: signature pass over the scene, Representation3DState updates per color-over-all |
| `bench-sig2.ts` | P12: split of the pass (signature, schema linear textures, key) |
| `bench-cuts.mjs` | P15: `targets.some` vs `Set` gathering at K=300/1500/3000 |
