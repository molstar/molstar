# Animation: `interpolate` nodes → frames

`generateStateTransition` (`helpers/animation.ts`) turns a snapshot's animation tree into `[MVS tree, durationMs]`
pairs, one per distinct frame; `assignStateTransition` (`load.ts`) converts and loads each ([load.md](load.md)).

## 1. Semantics

- **Steps**: `t = i × frame_time_ms`, `i = 0..ceil(duration / frame_time_ms)`; `duration` is the animation's
  `duration_ms` or the latest end of its transitions.
- **Groups**: transitions with the same `target_ref` and `property` (joined with `.`), sorted by `start_ms` then
  `duration_ms`. At `t` the active one is the last started (the first before any). A missing start value is the previous
  transition's end, else the property's value in the snapshot tree.
- **Target**: the node with that `ref`; `property` is in its `params`, or in `custom` when the first element is
  `'custom'`. A group whose node, `params` or `custom` is missing is skipped. The path is the first transition's
  `property`.
- **Each step starts from the snapshot tree**: a group that writes nothing leaves the snapshot's value.
- **A color with no start value** (no `start`, no previous end, nothing at the target): nothing is written.
- **Dedup**: a step equal to the previous one adds its duration to the previous frame. Noise (`Math.random`) makes every
  step distinct.

## 2. The update plan

The written paths are known before the first step: they are compiled once into an `ObjectMutator`
(`@molstar/core/util/object-mutator`) over the snapshot tree. A group is a slot (a leaf of a trie of keys); each step,
each group `set`s its value, `get` reads this step's value (set, else the snapshot's), and `commit` copies only the
objects on the paths of slots that changed since the previous commit, from it.

- **Dedup is identity**: `commit` returns the previous root when no slot changed. Slots whose path is missing in the
  snapshot tree ("dead") ignore writes.
- **Identity across frames**: a subtree no changed slot passes through is the same object in consecutive frames (and as
  in the snapshot tree); the conversion and loading memos rely on it ([load.md](load.md) §4).
- **Frames are immutable**: commits never write into an object an earlier commit returned; slot values are shared
  between frames.
- **Unset slots revert** to the snapshot's value (a key the snapshot lacks is removed).
- **Nested slots** (a slot inside another's value): the sets of a step are applied in order to a copy-on-write of the
  snapshot tree, dedup by `deepEqual`; no identity across frames.
- Slot equality is `deepEqual`, except distinct canonical values (§5), which are different.

## 3. Values

Values written into frames are fresh objects (vectors, matrices, maps); intermediate math uses module scratch.
Per-transition caches (palette function, rotation axis/angle, color map plan) are per group, by transition index.

## 4. Color maps

A transition of a color map (`group_colors`, a palette's `colors`) uses a `ColorMapPlan`: the keys in the order they are
set, each a constant or one of a few groups of keys sharing a start and end color. A step computes one color per group
and returns the map for those colors.

- Maps are kept by content per base map: equal maps from any transition of the animation are the same object, and
  canonical (§5).
- Plans: reused for equal raw start/end maps (`fromRaw`: same keys, order, strings), built in one pass when start and
  end cover the same keys (`fromSameKeys`), otherwise from decoded maps (`get`); equal plans are one object (`intern`).
- A step whose base map is not the plan's gets a new map over its own base.
- The caches (`ColorCaches`) belong to the animation's transitions and go with them after the expansion. Decoded colors
  are cached by string, module-wide and bounded.

## 5. Canonical values

`CanonicalValues` (`helpers/canonical-values.ts`): objects made once per content and never changed (the maps of §4). Two
distinct ones made for the same slot differ in content, so slot equality and `shareEqualParts` stop at identity
(`areDistinctCanonical`). Two from different slots may be equal: taking them as different can only miss reuse.
