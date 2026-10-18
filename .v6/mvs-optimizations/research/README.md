# MVS load optimizations: research

Goal: cut the time and memory of loading MVS states with animations — expanding the animations into frames and turning
each frame into a Mol\* state snapshot — to a fraction of today's, measured on real documents (motm1, exosome).

| File                                         | What it covers                                                                        |
| -------------------------------------------- | ------------------------------------------------------------------------------------- |
| [baseline.md](baseline.md)                   | fixtures, the bench, numbers, the profile by stage                                    |
| [animation.md](animation.md)                 | `helpers/animation.ts`: replacing mutative with a compiled update plan, allocations   |
| [molstar-conversion.md](molstar-conversion.md) | MVS → Molstar tree → state snapshot per frame: five steps, from local fixes to frames on demand |
| [specs.md](specs.md)                         | the `_spec` folders for `mvs-builder` and `mvs`, and the questions they must settle   |

## Findings

- exosome: 1.8 s and 3.3 GB allocated to load (Node, no canvas), 77 MB retained; motm1 0.74 s / 1.4 GB / 25 MB.
- 57 % of it rebuilds the whole Mol\* state tree of every frame through `StateBuilder`: UUID versions that are then
  overwritten by a params hash, a string hash per ref per frame, immutable.js maps, `toJSON`. A frame differs from its
  base snapshot in ~10 % of its transforms.
- 19 % is mutative (`produce`) in the expansion: proxies along every path every time step, then `deepEqual` to find
  the ~60 % of steps that are duplicates.
- 9 % is converting each frame's MVS tree: three full copies of the tree.
- Not in these numbers: a camera snapshot and a copy of the canvas props per frame when a canvas exists.

## Plan

| Step | Change                                                                                     | Where                    | Expected (exosome)                       |
| ---- | ------------------------------------------------------------------------------------------ | ------------------------ | ---------------------------------------- |
| 0    | goldens of today's output for the three fixtures; the bench (done)                          | `mvs/runtime/src/_test`  | —                                        |
| 1    | compiled update plan instead of mutative; value-level dedup (done)                          | `helpers/animation.ts`   | expansion 396 ms → tens of ms            |
| 2    | local waste: memoized ref hashes, no per-frame camera/canvas; cheap ids (done)               | `load*.ts`, `core/util`  | −25–35 %                                 |
| 3    | emit serialized transforms directly; reuse whole transforms from the previous frame (done)  | `load-generic.ts`        | most of the 57 %; retained 77 → few MB   |
| 4    | memoize conversion and loading by subtree identity, per unit (done)                         | `tree/`, `load-generic`  | per-frame cost ∝ changed units           |
| 5    | (open) frames computed at playback                                                          | `plugin/core` state      | animation load ≈ plan compilation        |

## Status

- **Step 0**: `load-animation.golden.test.ts` (in the unit run) and `load-animation.bench.test.ts`; how they work:
  `packages/mvs/runtime/_spec/performance.md`.
- **Step 1**: `ObjectMutator` (`@molstar/core/util/object-mutator`) replaces `produce`; dedup is root identity. exosome
  expansion 396 → 20 ms (795 → 26 MB allocated), motm1 140 → 10 ms; goldens unchanged. Precompiling easing and
  defaults was left: expansion is now 1.6 % of the load.
- **Ids** (from step 2): `UUID.create22` is a random 84-bit prefix per process plus a counter, ~70× faster (1M ids:
  515 → 7 ms); it was 8 % of the load. Still 22 URL-safe base64 characters, no longer random per call.
- Together: exosome `loadMVS` 1,803 → 1,243 ms, motm1 736 → 524 ms. Top of the profile now: `loadTreeVirtual` 66 %
  (`getTree` hashing 18 %, `getChildRef` 15 %, `StateBuilder.apply` 11 %), `convertMvsToMolstar` 12 %.

- **Steps 2–3**: `VirtualStateBuilder` (`load-generic.ts`) builds the serialized tree directly, checked against
  `StateBuilder`; ref hashes cached per load; `stringHash` one pass; frames skip camera, canvas and the `Entry` unless
  included; a frame reuses the transforms equal (as JSON) to the previous frame's or the snapshot's, and changed params
  share their equal parts. Goldens unchanged. exosome `loadMVS` ~800 ms (from 1,803), retained 19 MB (from 77); motm1
  ~270 ms (from 736), 4.9 MB (from 25). An API change after all: `UpdateTarget.update`/`commit` removed (ledger).
- **Step 4**: `convertMvsToMolstarMemo` memoizes conversion below the top nodes by MVS node; `SubtreeMemo` records the
  loading of each unit (a converted subtree right below a top node) and replays it while the unit object is the same:
  refs, transforms, serialized objects. Units are the granularity because the loading actions read their subtree; top
  nodes and units with `label_from_*` load every frame; a replay that does not fit falls back to a full load. Goldens
  unchanged, no fallback on the fixtures. exosome ~510 ms (61k units replayed, 20k loaded), motm1 ~200 ms.
- **Top nodes** (after step 4): `model`/`structure` collect annotations, tooltips and labels through a cache per unit
  object (`collectSubtree`); `nearestReprMap` is lazy; root children other than `download` (primitives, canvas,
  camera) are units; unit placeholders have a kind with no schema (the top conversion only copies them);
  `shareEqualParts` decides equality in its one walk, so a changed transform is compared and shared in one pass. Goldens
  unchanged. exosome ~510 → ~400 ms, motm1 ~200 → ~180 ms. A cache on every node made the non-memo path (features:
  frames with camera) 18 % slower: the caches are keyed by unit objects only.
- **Overall** (Node, no canvas): exosome 1,803 → ~400 ms (4.5×), 3.27 → 0.85 GB allocated, 77 → 19 MB kept; motm1
  736 → ~180 ms (4×), 25 → 4.9 MB.
- **Now on top** (exosome): the units loaded again (~7 of ~23 per frame: the components whose representation changed)
  — their actions ~60 ms and serializing their changed transforms ~100 ms; the top conversion ~45 ms; expansion ~20 ms.
  Next if needed: finer units (a representation's own inputs rather than its component's subtree), frames with
  `include_camera` (features) in the memo path by recording the camera context per unit, and a browser measurement.

Steps 1–3 were planned without an API change; step 3 removed `UpdateTarget.update`. Target after 1–3: exosome under 400 ms, motm1 under 150 ms, retained under 10 MB; the
bench gains a row per step. `_spec` files are written as each area is touched ([specs.md](specs.md)).
