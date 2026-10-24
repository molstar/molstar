# Slice 10: `prism/viewer/src` (v3)

The Prism viewer: API operations and their run path, the history's checkpoints, the sources (scene, style layers,
interaction, marks, sequence, state tree, clipping, structure move), the journal record, the Seq\* integration
(`sequence/`), plus the plugin's measurement sync and the callers of `Canvas3D.update(repr)`. HEAD `16277d822`, Node
26.12, Apple M4, headless (`createPrismViewerHeadless`, no rendering), N × `1cbs` and the 10-model `7qpd` (`multi.cif`).
The machine was loaded (load average 9–33 during the runs): wall times are ranges, quadratic claims rest on counters
and allocation. MVS is out of scope (`.v6/mvs-optimizations/`). P12 sits in frozen Canvas3D.

PRs since v2 that touched the slice: #2026 (P5 fixed, plugin side), #2032 (P14: `seekWhileViewing`), #2034 (B5:
`ChainInstance.span`), #2044 (P13/B6: frame step without clear and replay), #2037 (periodic trajectories:
`periodicByRef`), #2016 (bond providers in the source form), #2043 (sequences in the structure tree, track providers,
grouped structure tree, Model Archive PAE: `sources/sequence-scene.ts`, `sequence/model.ts` +618, `sequence/tracks.ts`,
`sources/scene.ts` +640, `molstar/secondary.ts`, `views/structure/sequence-inspector.tsx`), #2048 (headers only here).

## Hot paths now

| Path | When | Cost (headless, this run) |
| --- | --- | --- |
| `Canvas3D.update(repr)` → full precompile record (frozen) | per representation: move-gizmo drag, marks paint, plugin `SyncStructureRepresentation3DState`, model-archive pairwise | unchanged since v2: ~6.5 µs per render object in the scene per call, O(R_changed × R_scene) (P12) |
| structure tree's sequence rows (`refreshSequence`, #2043) | every sequence change: follow, track hide/show/add/remove, chains, style, letters, legend | one `row()` + `sequenceNode()` for **every** structure per bump: N=200 → 200–600 rows, 2.3–16.6 ms per op; N=50 → 0.5–1.7 ms (P17) |
| record append (`journal/record.ts` on `@molstar/ai` threads) | every recorded run, undo, turn delta | O(entries): 0.2 ms / 240 KB at 100 entries, 1.0 ms / 570 KB at 1000, 3.4–6 ms / 2.1–4.1 MB at 4000–8000 (P18) |
| single-structure op at N=200 (`representation.color`) | person/agent | 1.7–18 ms median (load 9–33; v2 2.9 ms at N=100 at load ~5); profile: checkpoint 7 %, hierarchy sync 7 %, record append ~6 % |
| frame step, sequence linked (7qpd, 20 steps) | playback, scrub | 0 plugin log lines (v2: 180), 0 record entries, 0 sequence syncs, 0 secondary-structure recomputes (default, `dssp` and `model` assignment); 1 view restart per step (P13 rest) |
| load at N=200 | per structure | 1 sequence row per load (no O(N) per load); Prism frames < 5 % of the load profile |
| group eye (Components, Measurements; #2043) | one click | one API run per row, in sequence: 13 rows 4–33 ms, 43 rows 18–48 ms (P19) |

## v2 → now

Counts: 11 v2 items open, partly fixed or fixed in v2's last pass (P5, P10, P12, P13, P14, P15, P16, B4, B5, B6, plus
P4 re-checked): fixed 5 (P5, P14, B5, B6; P4 stays fixed), partly 1 (P13), still open 5 (P12 frozen, P10, P15, P16,
B4), regressed 0, dropped 0.

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P12 | `Canvas3D.update(repr)` re-records every renderable's precompile signature | still open (frozen) | | [canvas3d.ts:2240](../../packages/graphics/src/canvas3d/canvas3d.ts#L2240) unchanged; no caller-side batching anywhere: [gizmos.ts:174](../../packages/prism/viewer/src/sources/structure-move/gizmos.ts#L174), [marks.ts:203](../../packages/prism/viewer/src/sources/marks.ts#L203) and [:213](../../packages/prism/viewer/src/sources/marks.ts#L213), [representation.ts:57](../../packages/plugin/core/src/behavior/static/representation.ts#L57)/63/69; one more caller found: [pairwise/ui.tsx:717](../../packages/extensions/model-archive/src/quality-assessment/pairwise/ui.tsx#L717) (a loop over every Representation3DState, moved in #2027) |
| P13 | Frame step restarts the link views | partly fixed | #2044 | no `highlight.clear`, no replay (`stop(false)`, `start(…, false)`, [link.ts:110](../../packages/prism/viewer/src/sequence/link.ts#L110)); still a new `linkViews`, `showMarks` (marks channels cleared and set again) and label provider per step: 20 `frameRebind` for 20 steps |
| P5 | (plugin) measurement sync walks the tree per change | fixed | #2026 | group ref from `cell.created`/`stateUpdated`, O(1) check per sync; reviewed: correct for create/remove in one transaction and for snapshot loads |
| P14 | `trajectory.frame` logs ~9 transform lines per seek | fixed | #2032 | `seekWhileViewing` with `doNotLogTiming` ([trajectory.ts:231](../../packages/prism/viewer/src/sources/trajectory.ts#L231)); 20 steps: 0 log events (v2: 180) |
| P15 | `findCuts`: linear `.some` per target, O(K²) | still open | | [cuts.ts:102](../../packages/prism/viewer/src/sources/clipping/cuts.ts#L102) unchanged (header only); v2's 8 ms at K=1500 stands |
| P16 | bulk runs' `did` one line per structure, kept uncut | still open | | [select.ts:351](../../packages/prism/viewer/src/api/select.ts#L351), [record.ts:520](../../packages/prism/viewer/src/journal/record.ts#L520) unchanged; now also costs the record append (P18) per entry size |
| P10 | `StateTreeSource` rebuilds on each unwatched read | still open | | [source.ts:318](../../packages/prism/viewer/src/features/state-tree/source.ts#L318) unchanged; no hot caller |
| B4 | outside changes never invalidate checkpoints while any animation plays | still open | | [checkpoint.ts:205](../../packages/prism/viewer/src/sources/checkpoint.ts#L205), `isViewing` ([trajectory.ts:240](../../packages/prism/viewer/src/sources/trajectory.ts#L240)) unchanged: camera spin still exempts every outside commit |
| B5 | `wholeChains` drops trailing hydrogens | fixed | #2034 | `chains.set(c, c.span)` ([link-views.ts:109](../../packages/prism/viewer/src/sequence/molstar/link-views.ts#L109)), span from the unit's chain atom segment |
| B6 | frame step clears the 3D hover and records `select.set` | fixed | #2044 | 20 steps: 0 record entries, 0 `highlight.clear` runs |
| P4 | (plugin) StructureInfo per created object | fixed (re-checked) | #2026 | `rootObjects` 0.9 % self of a load profile at N=200 |

New bugs or costs from fixes since v2: none from #2026/#2032/#2034/#2044. #2043 brought P17 (every structure's
sequence row per sequence change), P19 (group eye one run per row) and B7 (a group row's eye goes stale). #2043's
per-frame path is clean (no assignment flips, no recompute per step, measured).

## Open performance findings (ranked)

| Rank | Id | Status | Title | Where | Kind | Estimate | Effort | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P12 | still open | `Canvas3D.update(repr)` re-records every renderable; callers loop per representation | [canvas3d.ts:2240](../../packages/graphics/src/canvas3d/canvas3d.ts#L2240) | time | v2 (estimated from the signature pass): drag at N=100 ~12 → <0.5 ms per move; color over all O(N²): ~0.3 s at N=100, ~8 s at N=500 | S | yes (caller batching no) |
| 2 | P18 | new | Session record append copies and diffs the whole thread per entry: O(entries) per run, O(E²) per session | [host.ts:166](../../packages/ai/src/threads/host.ts#L166), [state.ts:51](../../packages/ai/src/threads/state.ts#L51), [state.ts:105](../../packages/ai/src/threads/state.ts#L105) | time+memory | measured: 1.0 ms / 570 KB per append at 1000 entries, 3.4–6 ms / 2.1–4.1 MB at 4000–8000 → ~0.2 ms / ~0.2 MB flat | M | no |
| 3 | P17 | new | Every sequence change rebuilds every structure's Sequence row | [scene.ts:1608](../../packages/prism/viewer/src/sources/scene.ts#L1608), [sequence-scene.ts:42](../../packages/prism/viewer/src/sources/sequence-scene.ts#L42) | time | measured: N=200, 200–600 rows, 2.3–16.6 ms per follow/track/style op → 1–2 rows (~0.1 ms) for per-structure changes | S | no |
| 4 | P15 | still open | `findCuts` O(K²) | [cuts.ts:102](../../packages/prism/viewer/src/sources/clipping/cuts.ts#L102) | time | v2 measured: 8 → 0.1 ms at K=1500 per state change while shown | S | no |
| 5 | P13 | partly fixed | Frame step still restarts the link views, marks channels and the label provider | [link.ts:106](../../packages/prism/viewer/src/sequence/link.ts#L106) | time | estimated: a view, 3 subscriptions, a label provider and 2 × marks channels per step fewer | S | no |
| 6 | P19 | new | A group's eye runs one operation per row, in sequence | [scene.ts:1226](../../packages/prism/viewer/src/sources/scene.ts#L1226) | time | measured: 43 rows 18–48 ms (13 rows 4–33) → one run | S | no |
| 7 | P16 | still open | bulk `did` one line per structure, uncut in the record | [select.ts:351](../../packages/prism/viewer/src/api/select.ts#L351) | memory | v2: ~10 KB per entry at N=200 → O(1) | S | no |
| 8 | P10 | still open | `StateTreeSource` rebuilds per unwatched read | [source.ts:318](../../packages/prism/viewer/src/features/state-tree/source.ts#L318) | time | plausible, no hot caller | S | no |

### P18: the session record's append is O(entries)

- **Mechanism**: `createSessionRecord` writes every recorded run as `host.append` on an `@molstar/ai` thread
  (thread made at [record.ts:384](../../packages/prism/viewer/src/journal/record.ts#L384)). The thread host's `run`
  ([host.ts:166](../../packages/ai/src/threads/host.ts#L166)) makes `draftOf(state)` (new `Map`s of every entry,
  stream, delta batch and receipt, a copy of `order`), runs the command, then `writesOf(state, draft)` compares every
  key of all four maps to find the 2–3 changed rows. The record is unbounded (append-only, no entry cap), so a session's
  k-th step costs O(k): a session of E steps is O(E²). The viewer's default rows are memory (`createMemoryRows`), so
  this is pure CPU and garbage on the main thread, after every run.
- **Evidence**: `bench-record.ts` (the real `createSessionRecord`, a stub API emitting `select.set` runs; 100 appends
  timed at each size): 202 µs / 239 KB at 100 entries, 308 µs / 346 KB at 500, 1031 µs / 569 KB at 1000, 3438 µs /
  1061 KB at 2000, 4155–6044 µs / 2127–2158 KB at 4000, 6052 µs / 4063 KB at 8000. CPU profile at 4000: `draftOf`
  35.9 % self, `writesOf`'s map loop 35.1 %, GC 7.5 %. In the viewer at N=200 (~250 entries) the thread code is ~6 %
  of a `representation.color` profile.
- **Fix** (in `packages/ai/src/threads`, the owner; the record needs no change): let commands write through a draft that
  records the keys it sets (copy-on-write maps over the loaded state, or a `touched` set per table), and make
  `writesOf` iterate only those; `order` appended, not copied. Keep the conflict retry (it reloads). Test: append cost
  and allocation flat from 100 to 8000 entries (a heavy-lane bound, or the scaling harness).
- **Risk**: the thread conformance suite (`packages/ai/src/testing/threads/conformance.ts`) covers the semantics.

### P17: every sequence change rebuilds every structure's Sequence row

- **Mechanism**: `createSceneSequence` bumps one `version` for any of the model's `tracks`, `edits`, `providers`,
  `theme`, `style`, `legend`, the followed structure, the band shown
  ([sequence-scene.ts:42](../../packages/prism/viewer/src/sources/sequence-scene.ts#L42)). The scene's
  `refreshSequence` ([scene.ts:1608](../../packages/prism/viewer/src/sources/scene.ts#L1608)) then calls `row()`
  (`sequenceChains`, `tracksOf` with `applying` sorting the providers, `addableTracks`, `describeFilter`) and
  `sequenceNode()` (menu closure, children, inspector update) for every structure in `structureRows`. Most changes
  touch one or two structures: an edit names one structure ref, a follow changes the old and new shown one, the band's
  shown flag the shown one. "Bulk operations are one change" is kept (one pass a burst), but the pass is O(N) for an
  O(1) change.
- **Evidence**: `bench-count.ts` (counters injected at bundle time by `build-inst.mjs`, sources untouched), per op:
  N=50: `sequence.follow` 2 bumps / 100 rows / 1.6 ms, `sequence.track.hide` on another structure 1 / 50 / 0.7 ms, on
  the shown one 3 / 150 / 1.3 ms, `sequence.style` 1 / 50 / 0.5 ms. N=200: 2 / 400 / 5.6 ms, 1 / 200 / 6.3 ms,
  3 / 600 / 16.6 ms, 1 / 200 / 2.3 ms. ~11–30 µs a row (loaded machine). Loads are fine: 1 row per load at N=200,
  `color`/`select`/`visibility` 0–0.5 rows.
- **Fix**: bump with a scope: `edits` diffed by structure key (the model already keeps them by ref), follow and band as
  {old, new} shown refs, `tracks` the shown ref; `providers`, `theme`, `style`, `legend` stay global (every row shows
  them). `refreshSequence` visits only the dirty refs (or all when global). Test: a `test:scaling` case, track hide at
  n and 4n structures, rows visited constant.

### P19: a group's eye is one run per row

- `groupEye` ([scene.ts:1204](../../packages/prism/viewer/src/sources/scene.ts#L1204)) runs `component.visibility` or
  `measure.visibility` once per row, awaited in sequence inside one `api.group` (one history step, K runs and K
  visibility updates). `bench-group-c.ts`: 13 components 4–33 ms, 43 components 18–48 ms, 13 / 43 runs. Each
  `measure.visibility` also re-reads the measurement list (`oneMeasurement` → `readMeasurements`), so Measurements is
  O(K²) reads (not measured headless: measurement labels need `canvas`). Fix: one operation over a list (a
  `components`/`measurements` param, `toggleVisibility([...])` once, `setSubtreeVisibility` per ref in one pass). Low:
  K is small for components.

### P13 (rest)

- The frame path ([link.ts:106](../../packages/prism/viewer/src/sequence/link.ts#L106)) still does `stop(false)` +
  `start(structure, view, false)`: `stopViews` tears down `linkViews` (three channel subscriptions), `showMarks` (clears
  every marks highlight channel and sets it again, `seq.highlight.clear/set` per mark) and the hover label provider, and
  makes them again for the new `SequenceView3D`. Fix as in v2: a `rebind(structure)` on the views, or keep `showMarks`
  and the label across a rebind (they read the link, which already follows the new Structure).

### P12 (frozen)

- Unchanged; v2's evidence stands. The non-frozen interim is still open: the drag could `setState` on every drawn
  representation then call `canvas3d.update()` once; `paintMarks` likewise once per paint; the plugin behavior once a
  task. The model-archive pairwise loop ([ui.tsx:717](../../packages/extensions/model-archive/src/quality-assessment/pairwise/ui.tsx#L717))
  is the same pattern over every Representation3DState.

## Bugs

| Id | Severity | Status | Title | Where | Failure scenario |
| --- | --- | --- | --- | --- | --- |
| B4 | medium | still open | While any animation (camera spin) or snapshot plays, outside changes never invalidate checkpoints | [checkpoint.ts:205](../../packages/prism/viewer/src/sources/checkpoint.ts#L205) | Camera spin plays; a script commits `plugin.build()…commit()`; Undo restores the pre-spin checkpoint and the edit is lost without the notice |
| B7 | low | new | A group row's eye (Components) is drawn from `hasEye` but nothing tells the row when it changes | [scene-tree.tsx:752](../../packages/prism/viewer/src/views/structure/scene-tree.tsx#L752), [scene.ts:331](../../packages/prism/viewer/src/sources/scene.ts#L331) | A structure with one shown component; hide its contents in Entities: the component becomes "emptied" (`hasEye` false), the group's `hasEye` turns false, but its label, description ("1 component") and hidden flag stay, so the row does not render again: its eye stays, and a click does nothing (`toggleVisibility` checks `hasEye`). The reverse: an eyeless group keeps no eye after its component returns |

- **B4**: as in v2 (code unchanged): exemption by time (`isViewing` = viewing counter or `isAnimating`; snapshot
  `isPlaying`), not by origin. Fix as in v2: exempt only the frame-changing animations' own updates.
- **B7**: `Node.hasEye` is a getter over the children's `hasEye` and `emptied`
  ([scene.ts:331](../../packages/prism/viewer/src/sources/scene.ts#L331)); the tree row's watch sources are label,
  description, hidden, picked and place ([scene-tree.tsx:752](../../packages/prism/viewer/src/views/structure/scene-tree.tsx#L752),
  [:1061](../../packages/prism/viewer/src/views/structure/scene-tree.tsx#L1061)), and the inspector reads `hasEye` once
  when it is built (`structure-tab.tsx`, `node.hasEye === false ? null : <VisibilityToggle>`). `emptiedNode` sets
  `emptied = true` ([scene.ts:935](../../packages/prism/viewer/src/sources/scene.ts#L935)) on the component's own node
  (same id, same count in the group's description), and `cg.set('Components', …, allHidden(comps))`
  ([scene.ts:1325](../../packages/prism/viewer/src/sources/scene.ts#L1325)) leaves the group's hidden flag as it was
  when the emptied row was the only one (`allHidden` of no eyed rows is false). Confirmed from code, not run in a DOM.
  Fix: make `hasEye` a state on `Node` (set where `emptied` and children change; groups derive it) and add it to the
  row's sources and the inspector. Test: the scene-tree DOM test with a component emptied and back.

## Already good, re-checked

- #2043's frame path: `bench-frames.ts` (7qpd, 20 steps, linked), default, `dssp` and `model` secondary assignments:
  0 `syncNow`, 0 `molstar:secondary_structure` recomputes, 0 sync DSSP attaches, 0 log events, 0 record entries.
- `sequenceChains` is cached per Structure and per hierarchy + chain set; O(units) per new Structure, a walk only of
  merged units.
- `componentKinds` (8 keyword queries) is memoized per Structure and runs only when the Components menu is read.
- `secondary.ts` `drawnAssignment`: reads the 3D view's computed property; attaches (synchronously) only when the 3D
  view has none, the same computation the cartoon would make.
- `sequence.ts`'s new `object.updated` subscriber: one `instanceof` and root check per update; the assignment string
  compare keeps the model quiet.
- The scene's `cell.stateUpdated` group refresh: O(rows of one group) per event; a bulk hide of K rows is K² row reads
  of `hasEye`/`isHidden`, negligible at the group sizes seen (≤ 50).
- The sequence model's `syncNow` (#2043): a handful of `findRef` walks over a document of tens to a few hundred nodes,
  once per change, serialized by `syncing`.
- #2026's measurement-group tracking: correct for a group removed and re-created in one transaction and for a state
  replaced by a snapshot (removed → stale → one search; created → tracked).

## Not audited

- The views' DOM cost in a browser (sequence inspector, grouped tree): jsdom and browser not run; the waste and
  subscription budgets (`waste.test.tsx`, `subscriptions*.test.tsx`) were not re-run.
- Measurements group eye at K > 0 headless: measurement labels need the `canvas` module, absent here.
- Model Archive PAE operations and card (`packages/extensions/model-archive`): extension slice; only its pLDDT track
  provider (`applies` is a property check) was read.
- P12 in a browser (no native GL here), as in v2.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/10/`:
`node build.mjs X.ts X.mjs` (esbuild, `molstar-src`) or `node build-inst.mjs X.ts X.mjs` (same, with counters injected
into `sources/scene.ts`, `sources/sequence-scene.ts`, `sequence/model.ts`, `sequence/link.ts`,
`sequence/molstar/secondary.ts` at bundle time); `node prof.mjs X.cpuprofile X.mjs.map` for profiles.

| Script | What |
| --- | --- |
| `bench-load.ts` | N × 1cbs loads with the sequence in; `TREE=0` without Sequence nodes; `PROF=load|color|follow` profiles |
| `bench-count.ts` | P17: sequence bumps, rows and `refreshSequence` ms per load and per op (`build-inst.mjs`) |
| `bench-record.ts` | P18: record append µs and KB allocated at 100–8000 entries |
| `bench-group-c.ts` | P19: the Components group's eye over K components |
| `bench-frames.ts` | P13, P14, B6, #2043 frame path: 20 steps on 7qpd, linked; `SS=dssp|model` |
