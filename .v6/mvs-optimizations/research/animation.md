# Expanding animations: `helpers/animation.ts`

`generateStateTransition(snapshot)` turns the `interpolate` nodes of a snapshot's animation into a list of MVS trees,
one per distinct frame, each with a duration. 19 % of exosome's load, 795 MB allocated for 2,886 frames.

## What it does now, per time step

For `i = 0..N` (`N = ceil(duration / frame_time_ms)`):

1. `produce(snapshot.root, draft => …)` (mutative): for each transition group (same `target_ref` + `property`),
   - `select(draft, nodePath)` walks `children[i]…` through draft proxies (a proxy per object on the path, each frame),
   - picks the active transition of the group by `start_ms` (linear scan from the end),
   - reads the start value from the **draft** when `start` is not given (`select(target, property)` again),
   - computes the value and `assign`s it through the proxies.
2. `deepEqual(root, prevRoot)`: when equal, the previous frame's duration grows instead.

## Costs, in order

| Cost                                                                                             | Fix                                                                         |
| ------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------- |
| mutative: a proxy and a shallow copy for every object on every path, every step; `finalize` walk | compiled update plan (below): plain copies along the dirty paths only      |
| node paths resolved through drafts (`select` 5 % self of the whole load)                         | resolve once per snapshot into the plan                                     |
| every time step builds a tree, duplicates are found afterwards by `deepEqual`                    | compare the **values** (slots) with the previous step; build nothing when equal |
| start values read from the tree every step                                                       | read once from the source tree, per transition                              |
| active-transition search, easing lookup (`EasingFnMap[…] ?? …`), `?? 0` defaults every step      | precompile each group: sorted start times, easing functions, defaults      |
| `Mat4()`, `Mat3()`, `Vec3()`, `Array.from({length}).fill(0.1)`, color strings per step            | allocate only for values that end up in a frame; reuse scratch for the rest |
| `getTransitionKey` builds strings; `Math.max(...transitions.map(…))`                             | once per snapshot (already), fine                                            |

## The update plan (the "ObjectMutator")

The paths an animation writes are known before the first step: `nodeMap[target_ref]` + `params`/`custom` + `property`.
So compile them once into a trie, and per step only write values into slots.

```ts
// once per snapshot
const plan = TreeUpdatePlan.create(snapshot.root);
plan.pushPrefix(...nodePath, 'params'); // e.g. 'children', 0, 'children', 3, 'params'
const opacity = plan.slot('opacity'); // slot = leaf index in the trie
const emissive = plan.slot(['custom', 'molstar_representation_params', 'emissive']); // relative to the prefix stack
plan.popPrefix();

// per time step
const m = plan.begin(); // values of the previous commit are kept
m.set(opacity, 0.42); // marks the slot dirty only if !== (or !equals for arrays) the previous value
const root = m.commit(); // the previous commit's root when nothing is dirty
```

- **Trie in flat arrays**: node `k` has `parent[k]`, `key[k]` (string or index), `firstChild/nextSibling` or a
  precomputed pre-order; leaves are slots. Per step: a `dirty: Uint8Array`, set on the slot and propagated to the root
  through `parent[]` (stop at an already dirty node).
- **Commit = path copy against the previous commit, not the source**: walk dirty trie nodes in pre-order, shallow-copy
  the object at that node (`Array.isArray ? slice() : { ...obj }`, or a per-node copier), write the child copies and
  the slot values. Everything not dirty keeps its identity from the previous frame. Two consequences:
  - frame dedup is `root === prevRoot`, `deepEqual` goes away;
  - unchanged subtrees keep their identity across frames, which the conversion can memoize on
    ([molstar-conversion.md](molstar-conversion.md) step 4).
- **Equality per slot kind**: numbers `===` (with `discrete` rounding applied before); vec3 / matrices: element-wise;
  colors: string `===`; color dicts: per key. Noise defeats it on purpose (every noisy step is a new frame).
- **Prefix API**: `pushPrefix`/`popPrefix` as proposed is for building the plan; per-step code holds slot numbers and
  never builds a path or a string.
- **Overlapping paths** (`color` and `color.5` in two groups; `matrix` and `matrix.3`): a slot inside another slot's
  value. Detect when compiling; write the outer slot first, then copy-on-write into the value just written (the trie
  node of the outer slot becomes an inner node for that step). Rare; a test for it.
- **Where it lives**: `@molstar/core/util/object-mutator` (done, `ObjectMutator`). Generic (nothing of MVS in
  it), hence in `core`. `@molstar/core/util/produce` stays for its other callers (`param-mapping`, the state builder).

## The rest of `animation.ts`

- Compile each group once: `targetSlot`, `kind`, transitions sorted by start with their `previous` resolved, the
  easing functions, `frequency`, `alternate`, `noise`, the decoded start/end (colors and rotations are already cached
  in `InterpolationCacheEntry`; make that the compiled transition rather than a `Map<transition, entry>` lookup per
  step).
- Values written into frames must be fresh (frames keep them): one array per changed vec3/matrix value per frame, no
  more. Scratch `Vec3`/`Quat`/`Mat4` for everything intermediate (`interpolateVec3`, `processTransformMatrix` already
  partly do).
- `interpolateScalars` for arrays: `new Array(n)` filled by the loop; the `.fill(0.1)` double-elements trick only
  matters for hot numeric arrays, not for values that go to JSON trees.
- Groups outside their time window write a constant (the start value before the first transition starts — which can
  differ from the node's own value — and the value at clamped `t = 1` after the last ends). Most groups
  of the long motm1/exosome animations are idle most of the time: compute the constant once and let the slot equality
  make the step free.
- `ctx.update` stays per step (it is cheap when `!ctx.shouldUpdate`).

## Semantics to keep (tests before the rewrite)

Golden output for the three fixtures: per snapshot, the list of `[frame tree, duration]` from the current
implementation, compared with the new one (`noise_magnitude` is random: seed `Math.random` in the test or skip noisy
groups). Plus unit tests for: `start_ms` ordering and the "previous end" start value, `frequency` + `alternate_direction`,
`discrete`, scalar ↔ array broadcasting, color dicts with base colors, palette colors, `transform_matrix` with pivot,
`custom.*` properties, a missing `target_ref` (ignored), overlapping paths.

## Expected

exosome 396 ms / 795 MB → a few tens of ms and a few MB (allocation ≈ the frames' own new values and path copies). The
larger effect is downstream: duplicate steps are never built, and unchanged subtrees keep their identity.
