# Animation: from `interpolate` nodes to frames

`generateStateTransition(ctx, snapshot, …)` (`helpers/animation.ts`) turns the animation tree of a snapshot into a list
of `[MVS tree, durationMs]`, one per distinct frame. `assignStateTransition` (`load.ts`) then converts and loads each
tree like a snapshot.

## 1. Semantics

- **Time steps**: `t = i × frame_time_ms` for `i = 0..ceil(duration / frame_time_ms)`; `duration` is the animation's
  `duration_ms` or the latest end of its transitions.
- **Groups**: transitions with the same `target_ref` and `property` (joined with `.`) form a group, sorted by `start_ms`
  then `duration_ms`. At time `t` the active transition is the last one started (the first before any has). A start
  value not given is the previous transition's end, else the property's value in the snapshot tree.
- **Target**: the node with that `ref` in the snapshot tree; `property` is in its `params`, or in its `custom` when the
  first element is `'custom'`. A group whose node, `params` or `custom` is missing is skipped. The group's path is the
  first transition's `property`.
- **Every frame starts from the snapshot tree**: a group that writes nothing at `t` leaves the snapshot's value.
- **Dedup**: a step whose tree equals the previous step's adds its duration to the previous frame. With noise
  (`Math.random`) every step is a frame of its own.

## 2. The update plan

The paths an animation writes are known before the first step, so they are compiled once per snapshot into an
`ObjectMutator` (`@molstar/core/util/object-mutator`) over the snapshot tree:

- the plan is a trie of keys (`pushPrefix` / `slot` / `popPrefix`); a group is a slot;
- each step, each group computes its value and `set`s its slot; `get` reads the value of this step (set, else the
  snapshot's) for start values and base colors;
- `commit` compares every slot with the previous commit (`deepEqual`) and copies only the objects on the paths of
  changed slots, from the previous commit. Nothing changed → the previous root itself.

Why: mutative's `produce` made a proxy per object on every path every step, then a `deepEqual` of whole trees found the
duplicates; 19 % of a load (exosome). The plan costs a few slot comparisons per step and builds nothing for a duplicate
step.

Rules that follow:

- **Dedup is identity**: `root === prevRoot`. Equal slot values ⟺ equal trees, except for "dead" slots (an object on the
  path is missing in the snapshot tree: writes are ignored, as `assign` did), which are left out of the comparison.
- **Identity across frames**: a subtree no changed slot passes through is the same object in consecutive frames (and in
  the snapshot tree until something under it changes). The conversion may memoize on it
  (`.v6/mvs-optimizations/research/molstar-conversion.md`, step 4).
- **Frames are immutable**: commits never write into an object a previous commit returned; slot values (arrays,
  matrices, color maps) are shared between frames and must not be mutated downstream.
- **Nested slots** (a slot inside another's value, e.g. `color` and `color.5` in two groups): the plan is not disjoint
  and the mutator falls back to applying the sets of a step in order to a copy-on-write of the snapshot tree, with
  `deepEqual` dedup against the previous frame — the old semantics exactly, without identity across frames.
- **Unset slots revert**: a slot not set in a step takes the snapshot's value; a key the snapshot does not have is
  removed again.

## 3. Values

Values computed per step go into frames, so each changed vector, matrix or color map is a fresh object; intermediate
math uses module scratch (`TransformState`, `RotationState`, …). Per transition caches (decoded colors, palette
function, rotation axis/angle) are per group, by transition index.

**A color with no start value** (no `start`, no previous end, the property missing — e.g. a key the target's color map
does not have): nothing to interpolate, the property is left as it is. Split coloring had made this an error ("Missing
start value for color interpolation"), which broke documents that relied on it.

## 4. Color maps

A color map transition (`group_colors`, a palette's `colors`: a dictionary of colors, thousands of keys in large
documents) is prepared once per content (`ColorMapPlan`): the keys in the order they are set, each a constant color or
one of a few groups of keys that go from the same start to the same end color. A step computes one color per group; the
map is looked up by those colors, then by the whole content (`ColorMapsByContent`, per base map), so equal maps from any
transition of a snapshot are the same object. Decoded colors are cached by string (a map repeats a few colors).

The plans, the raw maps they were made from and the maps by content are kept per animation (`ColorCaches`, held by its
transitions), so they are released after the expansion; only the maps the frames hold stay. A step whose base map is not
the plan's (a map slot inside another animated slot, where the base changes from step to step) gets a new map over its
own base, as every step did before plans.

Preparing a plan reads the raw start and end maps once: a transition whose maps equal those of an earlier one (same keys
in the same order, the same strings; documents repeat a few fades many times) takes that plan after one comparison
(`fromRaw`); maps over the same keys build it in one pass (`fromSameKeys`); plans with equal arrays are one object
(`intern`).

Why: a story flashing 10,596 line colors in 360 fades built a new map at every step of every fade and kept one in every
frame: 1.5 GB kept, most of a 39 s load. With one object per distinct map, slot comparison, frame dedup, hashing and
sharing stop at identity.

## 5. Canonical values

The maps of §4 are `CanonicalValues`: made once per content, never changed. Two different ones made for the same slot
have different content, so `areDistinctCanonical` lets the slot comparison (`ObjectMutator`'s `equals`),
`shareEqualParts` and `jsonEqual` answer without walking them. Two from different slots may be equal; taking them as
different costs only reuse (a new transform with the same version), never a wrong frame.
