# Load: MVS tree → state snapshot

`loadMVS` (`load.ts`) turns each snapshot, and each frame of its animation, into a `State.Snapshot` that is never
committed: the snapshot manager applies it later. Expanding the animation: [animation.md](animation.md).

## 1. Pipeline per snapshot and per frame

1. `convertMvsToMolstar` (`tree/molstar/conversion.ts`): defaults added, URIs resolved, MVS nodes split into Molstar
   nodes (`download`/`parse`, `trajectory`/`model`/`structure`), twins condensed. Three copies of the tree; a pure
   function of the MVS tree.
2. `loadTreeVirtual` (`load-generic.ts`): a DFS over the Molstar tree runs the loading action of each node kind
   (`MolstarLoadingActions`) and the loading extensions, which add transforms through `UpdateTarget`; the result is
   serialized.
3. For a snapshot: `molstarTreeToEntry` adds the camera (from `camera`/`focus` nodes), the canvas props, the duration,
   the entry. For a frame: only the state tree, plus camera and canvas with `include_camera` / `include_canvas`.

## 2. Refs

Transform refs are `!mvs:<stringHash(parent ref without !mvs: + ':' + transformer id[:repr type])>:<n>`
(`TargetManager`): the same node gets the same ref in every snapshot and frame, so applying a snapshot or a frame
updates cells instead of rebuilding them. The hashes are cached for the whole load (`refHashCache`): a frame computes
none. `stringHash` must keep its output (refs are in saved sessions); it is one pass over the code units.

## 3. The virtual builder and reuse

`UpdateTarget` builds on a `VirtualStateBuilder`, not a `StateBuilder`: a map of nodes with ordered children, no
immutable maps, nothing to commit. It reproduces the subset of `StateBuilder` the loader uses — a decorator transformer
is inserted between its parent and the parent's children; other transformers go under the parent's decorator root (the
end of a chain of single decorator children); tags and `dependsOn` are sorted unions; added and changed transforms get
the hash of their params as version (`useHashVersion`); the order is the pre-order with children in insertion order.
`_test/virtual-state-builder.test.ts` checks it against `StateBuilder`.

Why: building every frame through `StateBuilder` (immutable.js maps, a UUID version per transform overwritten by the
hash, `getTree`, `toJSON`) was 57 % of a load with animations.

A frame is serialized with `reuse = [previous frame, snapshot]`: a transform equal to the one with its ref there
(transformer, parent, tags, dependsOn, state, and params as `JSON.stringify` would write them, `jsonEqual`) **is that
object**, version included, and is not hashed. A transform that changed shares the parts of its params equal to the
previous ones (`shareEqualParts`). So a frame holds new objects only for what changed (~10 % of its transforms); the
rest are pointers into the previous frame and the snapshot.

Rules that follow:

- Serialized transforms in snapshots and frames are shared and must not be mutated (the plugin reads them; `fromJSON`
  makes new transforms around the same params objects, as with `StateTree.reuseTransformParams` before).
- Equality is `JSON.stringify` equality, not `deepEqual`: the version is a hash of `JSON.stringify(params)`, so key
  order, `undefined` and non-finite numbers matter. Transformers with `customSerialization` are never reused (their
  version hashes the raw params, which a serialized transform does not keep).
- The output is the same as before reuse: the goldens (`performance.md §2`) compare refs, params and versions.

## 4. Unchanged subtrees between frames

Frames of an animation differ in a few params; with the update plan ([animation.md](animation.md) §2) the MVS subtrees
nothing changed in are the same objects frame after frame. Two memos per snapshot use that:

- **Conversion** (`convertMvsToMolstarMemo`, `MolstarConversionMemo`): the top nodes (`root`, `download`, `parse` and
  their children) go through `convertMvsToMolstar` as before; a node below them converts to a copy of itself with
  defaults and resolved URIs whatever its parent or siblings (no rule, no condensing applies there), so it is memoized
  by MVS node. Its "units" are the converted subtrees right below the top nodes: a unit that is the same object as in an
  earlier frame has the same MVS source. `_test/conversion-memo.test.ts` checks it against `convertMvsToMolstar` on
  every frame of the fixtures.
- **Loading** (`SubtreeMemo`): loading a unit is recorded — the transforms it added (parent, transformer, ref, and after
  serialization what each serialized to), the ref counters it took, the MVS refs it mapped. When the same unit object
  comes again under the same parent target, the record is replayed: the same refs are taken, the transforms are placed
  with the builder's rules but serialize to the recorded objects, nothing is computed or compared.

Why units and not finer: the loading actions of a node read its subtree (children colors, transforms, clips), so a unit
reads only itself. What reads outside is excluded:

- the top nodes read their whole subtree (`structure`: labels, tooltips; `model`: annotations; `root`:
  `makeNearestReprMap`) — they are loaded every frame;
- `label_from_uri` / `label_from_source` read the nearest representation, which may be a sibling's — a unit containing
  one is loaded every frame (`isReusableSubtree`);
- a unit whose loading changes a transform it did not add (tags or state on its parent, a decorator inserted above
  others) is not replayable (`SubtreeRecord.reusable`);
- the camera and canvas are collected from the whole tree: no memo with `include_camera` / `include_canvas`;
- loading extensions must be `local` (read only their node, keep nothing in a context); one that is not turns the memo
  off.

A replay that does not fit (another ref counter, a replayed transform moved or changed) throws `SubtreeReuseMismatch`
and the tree is loaded again without the memo (`SubtreeMemo.stats.fallbacks`; the goldens require 0).
