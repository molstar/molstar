# Plan: a clip object dragged cheaply

## Status

- **Current**: steps 1 to 3 done; next: the yeast cell measured, then whether the structure representations want
  step 2. Benchmark case `clip/postsyn-drag` (PetWorld as the viewer adds its particles, culled on the GPU, a clip plane
  moved every frame; WebGPU, RTX): the representation's update a frame 1.59 -> 1.50 ms (step 1) -> 0.22 ms (step 2),
  the frame 2.72 -> 0.93 ms. PetWorld has few targets; the yeast cell's cost is per target, so it gains more.
  - **Step 1 done**: the merge check without `Object.keys` (a render object's values keep their keys: none
    merge-relevant is added after creation) and with the levels of detail compared as numbers in place.
  - **Step 2 done**: `BaseGeometry.updateClip` (the clip values alone) and one `Clip.getClip` per clip props object,
    shared by the visuals; the particle target representation sets the clip alone where it is the only change since
    its last full update (the same list, theme and other props) and keeps its defines (`isSameClipDefines`). Test:
    `repr/particles/representation/_test/target-clip.test.ts`.

## Question

Dragging a clip object (its preview, before the release) on a mesoscale model costs a long task each mousemove: on the
yeast cell about 74 ms (profiled, WebGPU), nearly all on the CPU. A drag only moves, turns or resizes a clip object:
what is drawn takes new clip uniforms, nothing else changes. Where does the time go, and how close to the uniforms'
cost can a preview get?

## Where the time goes

`ClipBinding.preview` (Prism, `sources/clipping/gizmo/binding.ts`) calls `repr.createOrUpdate({ ...props, clip })` on
every representation the cut is on. For the particle target representation (`repr/particles/representation/target/
representation.ts`) that is a full update of every target, per mousemove (profile, bottom-up, self / total):

- **Each target visual's `createOrUpdate`** (`target/common.ts`), 61 ms total: the props copied, the themes compared,
  the geometry props checked, nothing rebuilt; then `updateValues` and `updateRenderableState` over all its values
  (`ValueCell.updateIfChanged` 12.6 ms self, across every visual). All of it to set the clip values.
- **`setTransform` per visual**, 15 ms: returns early (the matrix is the same), once per target.
- **`mergeRenderObjectsByType`'s check that a merged group is still valid** (`isMergeFingerprintUnchanged`,
  `repr/merged.ts`), 41 ms total, 23 ms self, the largest: per member of every merged group an `Object.keys(values)`
  (a key count) and the levels of detail as a string (`lodOf`, 11.8 ms); both allocate (minor GC 9.5 ms).

Merged objects take their uniform cells from their first member (`createMergedValues`: `...first`): the members' clip
values set are what the merged objects draw with.

## Steps

1. **The merge check cheap** (every update of a merged representation, structures too): the `uLod` and `lodLevels`
   cells and their versions recorded as the other merge-relevant cells are, not a string; the `Object.keys` count
   dropped, or, if values gain keys after creation, found out without listing them. Expected: about half the time.
2. **A clip-only path in the particle target representation's `createOrUpdate`**: where only `clip` changed (the same
   particles and theme) and its defines stay (`dClipObjects`, the variant, the object count), each visual's clip values
   set and nothing else (no walk of the targets, no `setTransform`, no merge check). A drag keeps the defines; a change
   of them (a new shader variant, another merge) takes the full path. Then whether the structure representations (units,
   complex) want the same (far fewer visuals).
3. **Measured**: a benchmark case of a clip drag's preview (PetWorld, the yeast cell; `tests/src/benchmark/render.ts`):
   the time per mousemove and the frame, before and after each step.
4. **The GPU side**, if step 3 points there: a clip change invalidates the culling caches, so a drag culls again every
   frame (needed: the clip moved); barely visible in the profile.

Expected: a few ms a mousemove (a handful of value cell updates per visual) from about 70.

## Open questions

- Whether a merged object's members other than the first need their clip values during a drag (only the first's are
  read while merged; set all, so that an object leaving the group draws right).
