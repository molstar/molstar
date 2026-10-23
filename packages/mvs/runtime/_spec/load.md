# Load: MVS tree → state snapshot

Each snapshot, and each frame of its animation ([animation.md](animation.md)), becomes a `State.Snapshot` that is never
committed: the snapshot manager applies it later.

## 1. Pipeline

1. **Convert** (`convertMvsToMolstar`): defaults added, URIs resolved, MVS nodes split into Molstar nodes
   (`download`/`parse`, `trajectory`/`model`/`structure`), twins condensed (`download`, `parse`, `trajectory`, `model`).
   A pure function of the MVS tree.
2. **Load** (`loadTreeVirtual`): a DFS over the Molstar tree runs the loading action of each node kind
   (`MolstarLoadingActions`) and the loading extensions, which add transforms through `UpdateTarget`; then serialize.
3. **Entry** (`molstarTreeToEntry`): camera (from `camera`/`focus` nodes), canvas props, duration. A frame keeps only
   its state tree, plus camera and canvas with `include_camera` / `include_canvas`.

An animated snapshot and its frames go through the memos of §4; a snapshot without an animation does not.

## 2. Refs

`!mvs:<stringHash(parent ref without "!mvs:" + ":" + transformer id[":" + repr type])>:<n>` (`TargetManager`): the same
node gets the same ref in every snapshot and frame, so applying one updates cells instead of rebuilding them. `n` counts
refs with the same hash, in DFS order. Hashes are cached per load (`refHashCache`). `stringHash` must keep its output:
refs are in saved sessions.

## 3. The virtual builder

`UpdateTarget` builds on a `VirtualStateBuilder` (nodes with ordered children, nothing to commit) that reproduces the
part of `StateBuilder` the loader uses:

- a decorator transformer is inserted between its parent and the parent's children; any other goes under the parent's
  decorator root (the end of a chain of single decorator children);
- `tag` and `dependsOn` are sorted unions; `updateState` merges;
- serialized in pre-order, children in insertion order.

`_test/virtual-state-builder.test.ts` checks it against `StateBuilder`.

**Versions** (`setParamsVersion`, `jsonHash`): a hash of the params' JSON, so equal params have equal versions and
applying a frame updates only what changed. Params of up to 256 JSON entries hash their JSON text (= `hashMurmur128o`);
larger ones hash from the hashes of their parts, cached per object for the load (`jsonHashCache`), so parts shared
between frames are hashed once. Which way is decided by the JSON content alone.

**Reuse**: a frame is serialized with `reuse = [previous frame, snapshot]`. A transform equal to the one with its ref
there (transformer, parent, tags, dependsOn, state, params as JSON) **is that object**, version included:

1. the previous frame's is compared and its equal parts shared in one walk (`shareEqualParts`);
2. otherwise the version is computed (on the shared params), and the snapshot's transform is a candidate only with the
   same version, confirmed by `jsonEqual`;
3. otherwise the transform is new, with the shared params.

Rules:

- Serialized transforms are shared between snapshots and frames: never mutate them.
- Equality is JSON equality (`jsonEqual`), not `deepEqual`: the version hashes JSON, so key order, `undefined` and
  non-finite numbers matter. Transformers with `customSerialization` are never reused (their version hashes the raw
  params, which a serialized transform does not keep).
- `shareEqualParts` takes distinct `CanonicalValues` ([animation.md](animation.md) §5) as different without walking
  them: it can only miss sharing.

## 4. Replay of unchanged subtrees

The update plan keeps unchanged MVS subtrees the same objects from frame to frame ([animation.md](animation.md) §2). Two
memos per animated snapshot (`createFrameMemos`), shared by the snapshot and its frames, use that.

**Conversion** (`convertMvsToMolstarMemo`, `MolstarConversionMemo`):

- Top nodes (`root`, `download`, `parse`, the children of `download` and `parse`) go through the rules. Below them, and
  for the root's children other than `download`, a node converts on its own (no rule or condensing applies), so it is
  memoized by MVS node, and its converted `params`/`custom` by the MVS objects.
- Each `download` under the root is converted on its own and memoized, unless it has a twin to condense with or does not
  convert to exactly one Molstar download; then the whole top part converts as one.
- The **units** are the converted subtrees memoized this way (downloads, the subtrees right below top nodes).

**Loading** (`SubtreeMemo`): loading each node inside a unit is recorded: its own inputs, the transforms its action and
the extensions added (and what each serialized to), the ref counters taken, the MVS refs mapped, the target returned,
its children's records. A record is found by position (child indices from the root; the same in every frame) and, under
the same parent target:

- the same node object → the whole subtree is replayed: refs taken again, transforms placed by the builder's rules,
  serialized as recorded;
- the same own inputs (`isOwnReplayable` kinds only) → its own transforms are replayed and its target restored; its
  children are visited one by one;
- otherwise it is loaded and recorded.

**Own inputs**: a loading action, and an extension, read only their node (kind, params, custom, ref), its children that
have no action (colors, opacity, clips, transforms, instances, primitives, tooltips, labels) and the children's kinds;
never children with an action. A new loading action must keep to this, or its kind be excluded below.

Excluded:

- `model` and `structure` (`DeepKinds`) read their whole subtree, and `model` writes the annotation ids its descendants
  read: replayed only within a whole unchanged subtree, never own-only.
- `label_from_uri` / `label_from_source` read the nearest representation, possibly a sibling's: never replayed
  (`isReusableNode`), nor a subtree containing one.
- With `include_camera` / `include_canvas`, the nodes that write the camera or canvas into the context (`focus`,
  `camera`, `canvas`, `transition`) load every frame (`isReusableNodeWithContext`).
- A node whose loading changes a transform it did not add, and the node that added it, are not replayable.
- The memo is used only with `local` extensions and a `resetContext`.

A replay that does not fit (a ref counter, a moved or changed replayed transform) throws `SubtreeReuseMismatch`; the
context is cleared (`resetContext`: `MolstarLoadingContext.reset`) and the tree loaded again without the memo.

## 5. What top nodes read from units

Top nodes load every frame; some read their whole subtree, through caches:

- `model` (`collectAnnotationReferences`) and `structure` (`collectAnnotationTooltips`, `collectInlineTooltips`,
  `collectInlineLabels`): `collectSubtree` takes the items of a unit (`context.units`) from a cache keyed by the unit
  object. An item depends only on its node, its parent (inside the unit) and annotation ids (a function of params), and
  converted subtrees are never changed. Without `context.units`, nothing is cached.
- `root`: `context.nearestReprMap` is computed on its first `get`.

## 6. Mol\*-specific `custom` keys

Mol\* reads some keys of a node's `custom` that MVS does not define, prefixed `molstar_`: on a representation
`molstar_representation_params` (merged into its type params), `molstar_color_theme_name` and
`molstar_use_default_coloring`; on the root `molstar_on_load_markdown_commands`; on a shape `molstar_mesh_params`; on a
structure `molstar_bond_provider` (`{ name, params? }`, `structureProps`): the structure's bond provider, a name the
plugin lacks falling back to the built-in bonds with a warning when the snapshot is applied (bond-providers plan §4.2).
On a structure also `molstar_periodic` (`periodicModelProps`): a simulation's periodic box,
`'whole' | 'box' | 'nojump' | 'off'` or `{ mode, center? }` with center `'box'`, `'polymers'` or
`{ expression, language? }`; the conversion carries it to the structure's shown `model` node (not a topology's), whose
`ModelFromTrajectory` takes it as `periodic` and `periodicCenter`. Absent: molecules whole (the plugin's default); an
unknown mode is skipped with a warning. On a `coordinates` node `molstar_boxes` (`{ url }`): a NAMD `.xsc`/`.xst` file
giving the frames their periodic boxes with origins, by step; `CoordinatesWithBoxes` fetches it behind the coordinates
and carries the node's ref, and the conversion resolves the `url` against the view's as `url` params are
(`resolveCustomUrls`). A key is the way in for a Mol\* feature before (or instead of) a schema field.
