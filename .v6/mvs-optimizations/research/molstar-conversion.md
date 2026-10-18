# MVS tree → Mol\* state snapshot, per frame

`assignStateTransition` (`load.ts`) runs the whole single-snapshot pipeline for every distinct frame:

```
MVS tree ─ convertMvsToMolstar ─▶ Molstar tree ─ loadTreeVirtual ─▶ StateBuilder ─ getTree ─▶ StateTree ─ toJSON ─▶ State.Snapshot
           addDefaults (copy)                     LoadingActions        apply × n     useHashVersion   (serialized)
           resolveUris                            TargetManager refs                   (murmur of
           convertTree (copy)                     tags, dependsOn                       all params)
           condenseTree (copy + canonical JSON)
                                                  + camera snapshot, canvas props, Entry (UUIDs)
then StateTree.reuseTransformParams(frame, base)  // share params objects with the base snapshot when versions match
```

66 % of exosome's load (convert 9 %, `loadTreeVirtual` 57 %), and the retained 77 MB. A frame differs from its base in
~10 % of its transforms (baseline.md).

## Step 1 — local waste, no design change

| Waste                                                                                                                       | Fix                                                                                                            |
| --------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `TargetManager.getChildRef`: `stringHash(parentRef + ':' + suffix)` for every transform of every frame; 12 %                 | memoize `hashBase → hash` for the load (refs are deterministic: same parent and suffix give the same hash)    |
| `stringHash`: `input.split('').reverse().join('')` then two passes                                                         | one pass reading `charCodeAt` forwards and backwards; same output (refs end up in sessions: keep them stable) |
| `StateTransform.create` makes a `UUID.create22()` version that `useHashVersion` overwrites; `UUID.create22` is 8 % self      | a version option on create (or lazily), for builders that hash                                                 |
| `reuseTransformParams` builds the base's `Map<ref, transform>` again for every frame                                         | build once per snapshot                                                                                        |
| camera snapshot and `modifyCanvasProps(canvas3d.props)` (a copy of the canvas props) per frame, kept in the frame's entry    | compute only with `include_camera` / `include_canvas`; build no `Entry` for a frame                           |
| `validateTree(MolstarTreeSchema)` and other per-snapshot checks                                                              | already per snapshot, not per frame: keep it so                                                                |

Expected: 25–35 % of the load.

## Step 2 — emit serialized transforms directly

For a frame, nothing is committed: the `StateBuilder` (immutable.js maps of a transient tree), `getTree` and `toJSON`
exist only to produce `{ transforms: Serialized[] }`. A `VirtualTreeBuilder` behind the same `UpdateTarget` interface
(`apply`, `tag`, `dependsOn`, `deleteChildren`, `transformer`, `selector`) can push serialized transforms in pre-order
directly:

- `ref` from the memoized `TargetManager`; `parent` from the target; `transformer.id`; `params`; sorted `tags`;
  `dependsOn`; `state`.
- `version`: see step 3 — hashed only when the params differ from the reference transform.
- Pre-order is the builder's order of `apply` calls only if children are always applied after their parent and the
  extensions' `apply`s are placed the same way; to verify against `toJSON`'s pre-order (golden test: identical JSON
  for every frame of the three fixtures, base snapshots included).
- Use it for the base snapshot too (one code path); `replaceExisting` → the root's children are the new ones only, as
  now.

## Step 3 — share with the previous frame, not only with the base

Today a frame shares `params` objects with the base only (`reuseTransformParams`), and still holds 75 fresh serialized
transforms (with fresh `version` strings, tag arrays, …): ~27 KB a frame for exosome.

Instead, compare each emitted transform with the transform of the same `ref` in the previous frame (then the base):
params structurally equal (a deep compare of small JSON-like objects, cheaper than hashing them) and the same
`parent`/`tags`/`dependsOn`/`state` → reuse the **whole serialized transform object**, version included. Only
transforms that changed get a new object and a hash version. A frame becomes an array of ~75 pointers plus ~8 new
transforms: retained memory for exosome from 77 MB to a few MB, and no hashing for 90 % of transforms.

The format stays `PluginState.StateTransition` as it is; playback is unchanged (versions are equal where params are,
as now). Serialized sessions are unchanged (JSON expands shared objects, as it does today with shared params).

## Step 4 — recompute only what a frame changed (after the update plan)

With the update plan of [animation.md](animation.md), unchanged MVS subtrees keep their identity from frame to frame.
Then:

- **Conversion memo**: `addDefaults` per node is a pure function of the node: a `WeakMap<MVSNode, FullNode>` across
  the frames of a snapshot. `convertTree` rules read the parent (`parse` ← `download`, `structure` ← `parse`,
  `coordinates` ← `parse`): memo key includes the parent's identity for those kinds. `condenseTree` merges twins among
  siblings (`download`, `parse`, `trajectory`, `model`) — keyed by `canonicalJsonString(params)` — these nodes are
  never animated in practice, but the memo must be per sibling list, not per node.
- **Loading actions are not local per node**: `makeNearestReprMap` propagates a representation down to sibling
  subtrees (a `label_from_uri` can take the colors of a sibling's representation); `collectAnnotationReferences`,
  `collectInlineLabels`, `collectInlineTooltips` run at `model`/`structure` level over their subtree; `primitives`
  reads its whole subtree; the root action builds context for everything. The unit of reuse is therefore the
  **structure subtree** (or the root's direct children without a structure): re-run the loading actions of a unit
  only when its MVS subtree identity changed, else copy the unit's emitted transforms from the previous frame.
  Cross-unit reads (annotation sources by URI, primitives' `references`, `dependsOn` mapping) are refs, not params
  content: to confirm by a test that changes one unit and checks the others are bit-identical.
- Expected: the per-frame cost becomes proportional to the changed units (exosome: 66 animated targets in one snapshot
  but a handful of structures).

Do it only if steps 1–3 leave the conversion and the loading actions as the top of the profile.

## Step 5 (open) — frames on demand

The per-frame cost today is ~0.6 ms (exosome, Node); after steps 1–3 it should be well under 0.2 ms. That is far
inside a 16 ms frame, so frames could be **computed at playback** instead of at load: the transition keeps the compiled
plan and a cache of recent frames, load time for animations drops to the plan compilation. Costs:

- `PluginState.StateTransition` becomes (or gains) a provider, which is not JSON: session save must materialize the
  frames, or store the MVS source and regenerate on load (the MVS snapshot is already kept as `_transientData`).
- Frame boundaries (dedup) need the slot values of every step — cheap with the plan, without building trees.
- Scrubbing and `getStateTransitionFrameIndex` need durations upfront: from the same pass.

A plugin/core API change (`.v6/plans/breaking-v6-changes.md`); decide after steps 1–3 are measured, and with the
maintainers (it touches `PluginState` while the `PluginContext` rewrite is pending, Prism step 2).

## Equivalence tests

The steps above must not change the output. Before step 1, record goldens from the current code for the three fixtures:
per snapshot, the base `State.Snapshot` tree and every frame's `{ durationInMs, data }` as JSON (refs, transformers,
params, tags, dependsOn, state, and `version`: a hash of the params, so it must come out the same). Snapshot and entry
ids are UUIDs: left out. Kept as a test next to the bench, fixtures in `_test/fixtures`; the goldens are generated
from the current code once and committed (gzip if large), or computed from a frozen copy of the old pipeline.
