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
- Equality is `JSON.stringify` equality, not `deepEqual`: the version is a hash of the params' JSON, so key order,
  `undefined` and non-finite numbers matter. Transformers with `customSerialization` are never reused (their version
  hashes the raw params, which a serialized transform does not keep). Distinct `CanonicalValues` (§5 of
  [animation.md](animation.md)) compare as different without a walk.
- **Versions** (`setParamsVersion`, `jsonHash`): equal JSON, equal version, as `useHashVersion`. Small params (up to 256
  JSON entries) hash their JSON text, which is `hashMurmur128o`; larger ones hash from the hashes of their parts, kept
  per object for the load (`jsonHashCache`), so a part shared with an earlier frame (coordinate lists, a color map) is
  not walked again. Why: hashing the whole JSON of each changed transform was 39 % of a story whose primitives carry ~1
  MB of params (18 s of 46).
- The output is the same as before reuse: the goldens (`performance.md §2`) compare refs, params, versions, and the
  camera and canvas of frames.

## 4. Unchanged subtrees between frames

Frames of an animation differ in a few params; with the update plan ([animation.md](animation.md) §2) the MVS subtrees
nothing changed in are the same objects frame after frame. Two memos per snapshot use that:

- **Conversion** (`convertMvsToMolstarMemo`, `MolstarConversionMemo`): the top nodes (`root`, `download`, `parse`, the
  children of `download` and `parse`) go through `convertMvsToMolstar` as before, each unit a placeholder of a kind with
  no schema; a node below them, or a child of `root` other than `download`, converts to a copy of itself with defaults
  and resolved URIs whatever its parent or siblings (no rule, no condensing applies there), so it is memoized by MVS
  node, and its converted `params` and `custom` by the MVS objects: a node converted again because something below it
  changed keeps them. Its "units" are those subtrees (components under a structure; `primitives`, `canvas`, `camera`
  under the root), in the same order in every frame. `_test/conversion-memo.test.ts` checks it against
  `convertMvsToMolstar` on every frame of the fixtures.
- **Loading** (`SubtreeMemo`): loading each node of a unit is recorded — its own inputs (below), the transforms its
  action and the extensions added (parent, transformer, ref, and after serialization what each serialized to), the ref
  counters taken, the MVS refs mapped, the target its action returned, and the records of its children. A node finds its
  last record by position (the unit's index, then child indices) and, under the same parent target:
  - the same node object → the whole subtree is replayed: refs taken again, transforms placed with the builder's rules
    but serializing to the recorded objects, nothing computed or compared;
  - the same own inputs → only its own transforms are replayed, its target restored for its children, which are visited
    one by one (a representation whose emissive changed reloads alone, not its component);
  - otherwise it is loaded and recorded.

**Own inputs**: a loading action, and an extension, read only their node (kind, params, custom, ref) and the node's
children that have no action of their own (colors, opacity, clips, transforms, instances, primitives, tooltips, labels)
plus the children's kinds (`isPhantomComponent`); never children with an action (a component's representations). So a
node whose own inputs are the same objects adds the same transforms. A new loading action must keep to it, or its node
kind be excluded like `label_from_*` below.

What reads outside is excluded:

- the top nodes read their whole subtree (`structure`: labels, tooltips; `model`: annotations) — they are loaded every
  frame, from what §5 caches per unit;
- `label_from_uri` / `label_from_source` read the nearest representation, which may be a sibling's — loaded every frame
  (`isReusableNode`), and so is no subtree containing one replayed whole;
- a node whose loading changes a transform it did not add (tags or state on its parent, a decorator inserted above
  others) is not replayable, nor is the node that added it (`SubtreeRecord.reusable`);
- the camera and canvas are collected from the whole tree: no memo with `include_camera` / `include_canvas`;
- loading extensions must be `local` (read only their node, keep nothing in a context); one that is not turns the memo
  off.

A replay that does not fit (another ref counter, a replayed transform moved or changed) throws `SubtreeReuseMismatch`
and the tree is loaded again without the memo (`SubtreeMemo.stats.fallbacks`; the goldens require 0).

## 5. What the top nodes read from their units

The top nodes are loaded every frame, and some read their whole subtree. They read it through a cache or lazily:

- `model` (`collectAnnotationReferences`), `structure` (`collectAnnotationTooltips`, `collectInlineTooltips`,
  `collectInlineLabels`): `collectSubtree` walks the tree and takes the items of a unit (`context.units`) from a cache
  keyed by the unit object. An item depends only on its node, its parent (a tooltip's or label's component, inside the
  unit) and annotation ids (`stringHash` of the spec: a function of params), and a converted subtree is never changed,
  so the cached items stand. They are the same objects frame after frame: the `CustomModelProperties` and
  `CustomStructureProperties` params of two frames compare by identity below their wrappers.
- `root`: `context.nearestReprMap` is computed on its first `get` (a structure with inline labels, a `label_from_*`),
  not every frame.
- Outside the frame memo (`context.units` unset: snapshots, frames with camera or canvas) nothing is cached: every node
  would be new, and caching every node made those loads 18 % slower.

Why: with units replayed, these walks over every node were the largest per-frame cost left (exosome: ~140 ms of ~530).
