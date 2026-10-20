# 14: many separate structures in one scene

The other end from [13-large-scale](13-large-scale.md): up to 500 separate small and normal-sized entries (1.2k–11.5k
atoms each, 5.3k on average), each its own data, model and structure. Nothing is shared between structures, so the
scene holds hundreds of structures, thousands of state cells and representations, and the costs that matter are the
fixed ones per object: the plugin and state work per cell, structure and representation, and anything that walks all
of them for each one (O(objects²)).

## Bench

`tests/src/pipeline/many.ts` (run through `scripts/pipeline-bench.mjs many`; fixtures `many-fixtures.ts`: six local
entries cycled, each moved to its own grid cell). A real `PluginContext` with the default spec loads N files (raw data,
parse, the default hierarchy preset split in stages, the `auto` representation preset), then a theme change for all and
for one structure, 20 hovers and 10 clicks on one structure, overpaint of one and of all, visibility, removal.

- `--transaction stage` (the UI paths: `OpenFiles` and Prism's `load.files` wrap a load in one data transaction),
  `each` (one per file), `none` (plugin-API scripts: every commit notifies the managers).
- No WebGL in Node: `plugin.canvas3d` is a stub that does Canvas3D's work on representations (add, remove, mark every
  representation, update) and counts what reaches the frozen scene (`scene.update` over all render objects after a
  mark that changed something). The graphics context is bench.ts's CPU stub (impostors on). Canvas3D and gl internals
  are frozen and were not changed; their cost per object is not measured here.
- Node 25.6, Apple M4, `NODE_ENV=production` as the app builds; the desktop was in use (load 3–7), so runs were
  interleaved A/B, 2–4 rounds, medians; the allocation column is the stable one.

## What dominated (before)

N = 300, stage mode, `repr:auto` 25 s (84 ms per structure), allocating 27 GB; N = 100 → 300 the time per structure
×2.6 (none mode ×3.6, 114 GB allocated at N = 300). Geometry was 5–7 % of it. Ranked (share of `repr:auto`):

| #   | Cost                                                                                          | Share (stage / none) | Kind     |
| --- | --------------------------------------------------------------------------------------------- | -------------------- | -------- |
| 1   | StructureInfo: four `rootsOfType` walks of the whole state per created object                 | 40 % / 11 %          | O(N²)    |
| 2   | State update: `findDeletes` + `findUpdateRoots` walk both trees per multi-edit commit         | 31 % / 22 %          | O(N²)    |
| 3   | Measurement manager: six `findTagInSubtree` walks per state change                            | – / 38 %             | O(N²)    |
| 4   | Structure, volume, particle hierarchies rebuilt by a whole-state walk per state change        | – / 30 %             | O(N²)    |
| 5   | `updateFocusRepr` once per structure: focus behavior update, 3 walks + an empty commit        | 8 %                  | O(N²)    |
| 6   | SelectLoci marks the (empty) selection of a new representation on every representation       | 2–13 %               | O(N²)    |
| 7   | Geometry (cartoon trace mesh, impostors), bonds, parse                                         | 5–7 % → most         | per atom |

The none-mode rows only run outside transactions (scripts); in the UI paths they run once per action, O(N).

## Fixed in this branch

One commit each, identical results pinned by `_test/` tests recorded on the code before (14/F1–F8):

| Id  | Change                                                                                                              |
| --- | ------------------------------------------------------------------------------------------------------------------- |
| F1  | StructureInfo sets the max indices only after a root model or structure was created, updated or removed             |
| F2  | the measurement manager looks for its group once per change                                                         |
| F3  | the focus representation update walks once, commits only edits                                                     |
| F4  | state updates take deletes and roots from a builder's recorded edits (core; differential test against the old diff) |
| F5  | `updateFocusRepr` skips themes already set                                                                          |
| F6  | StructureInfo keeps its root models and structures from the events instead of walking                               |
| F7  | SelectLoci does not mark an empty selection on a new representation                                                 |
| F8  | the volume and particle hierarchies skip the walk without any volume or particle list                               |

Master meanwhile landed 10/P4 and 10/P5 (#2008), versions of F1 and F2 from the Prism audit; the merge keeps them:
the measurement manager is master's (it also skips the walk without a group), and StructureInfo has master's
semantics (a removed model lowers the max at once; a replaced or remapped object gets the last max and keeps its
index) computed over F6's tracked roots, so it walks for no object. The per-commit table below is from before the merge.

One difference besides speed (F3, F5): the focus behavior's empty data commit is gone, and it used to clear the data
state's undo history (an update without `canUndo`), toggle `isUpdating` and wait for queued updates. The presets'
own commits clear the history anyway.

### Per commit (N = 300, stage mode, ms / MB allocated, median of 2 interleaved rounds)

| Stage           | base         | F1           | F2           | F3           | F4          | F5          | F6          | F7          | F8          |
| --------------- | ------------ | ------------ | ------------ | ------------ | ----------- | ----------- | ----------- | ----------- | ----------- |
| model+structure | 2314 / 4622  | 1184 / 2356  | 1149 / 2348  | 1157 / 2342  | 1345 / 2355 | 1208 / 2344 | 203 / 200   | 188 / 202   | 163 / 186   |
| repr:auto       | 19167 / 27096 | 12360 / 14437 | 11790 / 14432 | 10818 / 11716 | 7163 / 3543 | 5969 / 2644 | 6420 / 2675 | 5469 / 2617 | 5358 / 2646 |
| overpaint:all   | 1946 / 3230  | 109 / 83     | 98 / 64      | 109 / 65     | 154 / 60    | 108 / 65    | 103 / 63    | 108 / 63    | 92 / 54     |
| theme:one       | 66 / 42      | 61 / 42      | 42 / 30      | 50 / 29      | 44 / 22     | 45 / 23     | 45 / 23     | 47 / 23     | 26 / 13     |
| remove:one      | 66 / 45      | 72 / 43      | 53 / 30      | 51 / 29      | 47 / 23     | 58 / 24     | 44 / 24     | 47 / 23     | 35 / 12     |

The desktop was busy: the time of `repr:auto` varies ±10–20 % between rounds (F6's is noise); the MB do not.
Measured before the review fixes of F4–F6 (pre-order index, overlapping `updateFocusRepr` calls, roots of the current
update); the final head against F8 (2 rounds): `repr:auto` 4607 → 4477 ms, 2606 → 2653 MB, model and structure
133 → 141 ms: the same within noise.

### Scaling (stage mode, before → after)

| N   | load to full scene¹ | `repr:auto` per structure | `repr:auto` allocated | overpaint of all  | theme / overpaint / removal of one |
| --- | ------------------- | ------------------------- | --------------------- | ----------------- | ---------------------------------- |
| 100 | 4.4 → 2.8 s         | 28.8 → 14.7 ms            | 3.8 → 1.0 GB          | 282 → 31 ms       | 20 / 20 / 20 → 16 / 9 / 8 ms        |
| 300 | 25.2 → 8.8 s        | 63.9 → 17.9 ms            | 27.1 → 2.6 GB         | 1946 → 92 ms      | 66 / 73 / 66 → 26 / 25 / 35 ms      |
| 500 | 74.3 → 16.3 s       | 120.9 → 20.8 ms           | 72.3 → 4.2 GB         | 6906 → 187 ms     | 122 / 147 / 435 → 45 / 46 / 46 ms   |

¹ load + model and structure + bonds + secondary structure + `repr:auto`. At N = 500, per event: hover 6.2 → 5.1 ms,
select 13.2 → 5.4 ms, deselect all 62 → 24 ms.

Without transactions (none mode, N = 200, plugin-API scripts): `repr:auto` 28.6 → 8.6 s (51 → 11 GB), model and
structure 3.4 → 0.9 s, overpaint of all 1001 → 65 ms, theme of one 31 → 13 ms.

## Not fixed (ranked)

1. **Hover, select and deselect visit every representation** (07/P3's other half). `LociHighlightManager` and
   `LociSelectManager` mark through the canvas without a representation, and Canvas3D's `applyMark` loops over all of
   them (`repr.mark` early-outs on `areRootsEquivalent`); a changed mark then updates every render object of the scene
   (`scene.update(void 0, true)`: 100k render-object visits for 20 hovers at N = 300). At N = 500 (1.6k
   representations) about 5 ms per hover and per click (bench, with its GC) before any GPU work. An index from root structure to
   representations in Canvas3D (frozen) and a scene update limited to the marked render objects would make both
   O(representations of the structure). Not measured in a browser.
2. **The structure hierarchy is rebuilt by a whole-state walk on every state change** (outside transactions): 12 ms
   at N = 500 of each single action (theme, overpaint, removal of one structure), plus one measurement walk (3 ms).
   An incremental hierarchy (refs updated from the cells an update touched) would make actions O(1) in N. Prism's own
   subscribers to state changes were not measured.
3. **Scripts that load without a transaction** pay every manager's walk per commit: none mode is still O(N²) (`repr:auto` at N = 200 8.6 s
   against 5.4 s at N = 300 in a transaction). Wrapping loads in `plugin.dataTransaction` (as `OpenFiles` and Prism do) is the cheap answer; worth a
   line in the plugin docs.
4. **Cartoon meshes reserve about 1.4× their vertices** (fill 0.70; sheets use 8 vertices per segment, the estimate
   `radialSegments`). `compactTrimmed` (08/P3, #2007) then copies 6.1k of 8.2k buffers (1.56 GB at N = 300, about
   0.35 s of `repr:auto`) and saves 690 MB retained, so it stays; a per-element estimate in `createPolymerTraceMesh`
   would avoid both. The item-count sizing of `CylindersBuilder`/`LinesBuilder` (08/P8) is negligible here (link
   cylinders 0.7M elements at N = 100, fill 0.53).
5. **Per-structure work that remains linear but heavy**: intra-unit bonds 6 ms per structure (`findBonds`, unit
   lookups), parse 4 ms, the cartoon trace mesh. Per atom, not per object; see 05/P7 and 08.
6. **13/P4** (a per-visual cache of the foreign-resource walk, branch `verify-p4-p7`) re-measured here: `theme:all` at
   N = 300 231 → 217 ms (−6 %, within noise), nothing else moved. With small visuals the walk is about a tenth of a
   recolor; it does not justify itself at this end.
7. **Bug found by the differential test** (pre-existing, core TODO): random builders with generated refs can leak a
   child added by a later builder into a child set of `state.tree` before that builder commits; the next whole-tree
   diff then fails in `findDeletes`. Not seen with explicit refs.

## Reproduce

```sh
node scripts/pipeline-bench.mjs many-fixtures --count 500
node scripts/pipeline-bench.mjs many --count 300 --transaction stage
node scripts/pipeline-bench.mjs many --count 300 --transaction stage --profile /tmp/p --only repr:auto
node scripts/pipeline-bench.mjs analyze /tmp/p 30
```
