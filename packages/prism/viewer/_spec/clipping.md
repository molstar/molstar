# Clipping

The `clipping` feature (`src/features/clipping`, in the default spec): cuts through representations (a plane, a sphere,
a box, a cylinder, a cone), made and edited in the Selection tab's Clipping section, from the selection menu's Clip… and
the scene tree's Clip this…, with handles in 3D, kept whole around chosen atoms, and drawn in the scene. Its plugin-side
code is in `src/sources/clipping` (the viewer's boundary: only sources and operations touch the plugin).

## 1. Cuts are read from the state

- **The state tree is the record.** A cut is a Mol\* clip object in representations' params (`type.params.clip` or
  `clip`: `RepresentationClip.get`/`set`); the feature keeps no store. A clip set from code (MVS, mesoscale, a state
  update) is listed too. Why: undo, sessions and snapshots hold cuts with nothing to sync.
- **One cut, several representations.** The feature tags its objects (`Clip` objects' hidden `tag`, `clip-<uuid>`, kept
  through edits): a tagged cut is the objects sharing its tag; an untagged one, those equal in every value
  (`RepresentationClip.sameCut`). Two cuts alike stay two when tagged.
- **Numbered in tree order** (`findCuts`): `Cut { index, object, targets, drawing?, kept? }`. The component manager's
  canvas-wide objects are not cuts; the cuts' drawings (§6) are never cut.
- **Per representation, not per atom** (D2): a cut applies to whole representations; atoms are left whole by §5.

## 2. The `clipping` area

- `list`; `add` (`shape`; `target` atoms, the selection by default, or `node` (volumes, shapes), or `everything`;
  `types`; `invert`; `keepWhole`); `update` (`shape`, `position`, `normal` or `rotation`, `scale`, `invert`: only what
  is given); `apply` (one or several representations on or off, one step; never the last off); `keep` / `release` (§5);
  `show` (§6); `remove`; `clear`. Each a person's step, each a command, in the `clipping` skill.
- **A new cut** is placed by its target's atoms (`placeCut`): a sphere or a box around a small target, with room; a
  window toward the viewer into a large one (one centered inside would show nothing); a plane through the center facing
  the camera; a cylinder a shaft along the view. With nothing selected, everything drawn.
- **Add again** returns the cut like it on those representations (shape, flip, place, direction, size within 0.5 Å and
  ~2.5°: `nearlyEqual`) instead of stacking another.
- **Another shape** keeps the place, the direction (a plane's normal is a cylinder's axis, a cone's opening), the size
  and the flip (`reshape`); a plane or cone has none, so the new one is sized from what it cuts.
- **Directions** are unit; one of length 0 is refused.

## 3. The section

- **Cuts · + Add**: a plane through what the selection is drawn as (narrowed by Through), or everything drawn; the new
  cut (or the one like it) becomes active. Keep the selection whole (§5) for new cuts.
- **One active cut**: the one whose editor is open (its name toggles it; Add and Clip this… make theirs), its handles
  always shown. Why: a cut drawn in the scene is never one that cannot be dragged.
- **A row**: the shape's icon and name, what it cuts (a count past three), drawn in the scene (the `scan` toggle),
  remove.
- **The editor**: the shape (a strip, as a representation's type), cut away the other side, then folds closed at first:
  Cuts through (x/N: by type, then each, those it is on first, searched past eight), Kept whole, Position and size, In
  the scene (while drawn). Forms send only the fields changed (a box's up alone would drop its twist), one update in
  flight.
- **Clip…** in the selection menu (`ui.selectionActions`) opens the section; **Clip this…** in a row's menu
  (`ui.nodeActions`: components, representations, volumes, particles) adds a plane through what the row draws, active,
  the selection left as it is.

## 4. Gizmos

- **Outside the state tree**: a `VisualSet` in the canvas; their clicks claimed (`PluginContext.clickOwners`: no focus,
  no select). Shown for the active cut while the section is shown (`onActive`).
- **Each handle says what a drag does** (`Handle`): move (the screen's plane), slide (along the axis; foreshortening
  corrected, capped at 4×), turn (0.5°/px about the screen's axes), size (whole: by the cursor's distance from the
  center over the grab's), faces (a box's face knob: its face, the opposite one stays), across (a box's four edges along
  an axis: the section about the center), radius (a cylinder's rims), angle (a cone's base rim). Outlines are handles; a
  sphere's ring faces the camera.
- **A drag runs live** on the representations (clip-only `createOrUpdate`: uniforms only, `repr/_test/clip-update`),
  coalesced, from the press point (CSS px); the release is one `clipping.update`, of the cut found by its object (its
  number may shift). Esc, a change of what it is on, or a failed commit cancels: the live representations go back to the
  state, the gizmos resync, the trackball's bindings are restored.
- **Bodies** are the clip volumes' walls (`graphics/geo/clip-volumes.ts`), drawn on both sides as fronts; handles glow.

## 5. Keep whole

- **Mol\*'s clipping layer** on a structure representation (the cell `representation.clipping` writes, tagged
  `clipping-controls`): group bit _i_ means "not clipped by clip object _i_". A cut keeps atoms by its object's index
  there. Only the first six (Mol\*'s groups).
- **Bits on the atoms' flags**: an atom's flags are those of the last layer covering it, so a bit is set or cleared on
  them (`withKept`, `withoutKept`), never by a layer over them.
- **Removal remaps** in the same step: a cut taken off a representation drops its bit and the later bits move down
  (`afterRemoval`); the component manager replaces its own objects in place (`replaceOwned`). Why: an index is what a
  group means.

## 6. Drawn in the scene

- **`ClipVolumesShape3D`** (plugin core, imported by the feature: registered on import, Prism only): from a
  representation that can be cut to a shape provider of its clip objects (by tag) as translucent walls and solid edges
  (`ShapeProvider.transparency`, per group); options for walls, edges, colors, opacity, width, extent (a plane's disk
  and a cone's reach relative to what they cut, 0.25–4×: in proportion when that changes). Built from its data (objects,
  bounds, options): the shape representation keeps its first `getShape`.
- **One drawing a cut**, deduplicated by its tag (`clip-volume-<cut tag>`; every cell also `clip-volumes`), under the
  cut's first representation; it moves when the cut leaves that one and goes with the cut. Saved with the state.
