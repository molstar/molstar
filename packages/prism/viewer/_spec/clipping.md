# Clipping

The `clipping` feature (`src/features/clipping`, in the default spec): cuts through representations (a plane, a sphere,
a box, a cylinder, a cone), made and edited in the Selection tab's Clipping section, from the selection menu's Clip… and
the scene tree's Clip this…, with handles in 3D, kept whole around chosen atoms, and drawn in the scene. Its plugin-side
code is in `src/sources/clipping` (the viewer's boundary: only sources and operations touch the plugin).

## 1. Cuts are read from the state

- **A sphere is an ellipsoid**: the shaders clip it so (its half axes the scale's halves, turned by its rotation), so it
  is drawn, sized and edited so everywhere (walls, edges, gizmo, form); round when placed. The type stays `sphere`
  (Mol\*'s, MolViewSpec's, saved states). Why: a scale set from code or MolViewSpec cut another volume than was drawn.
- **A cylinder is elliptic**: Mol\*'s clip takes its radii from the scale's x and z halves (its height y), and the cut
  is drawn and edited so; round when placed. The cone stays round (its angle alone).

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
  `types`; `invert`; `keepWhole`); `update` (`shape`, `position`, `normal` or `rotation`, `twist`, `scale`, `invert`:
  only what is given); `apply` (one or several representations on or off, one step; never the last off); `keep` /
  `release` (§5); `show` (§6); `remove`; `clear`. Each a person's step, each a command, in the `clipping` skill.
- **A cut by its number or its tag** (`cut`, `tag`): the views and the gizmos send the tag. Why: runs wait their turn,
  so a number read at a click can name another cut by the time the run starts (a double click's second remove).
- **Params checked**: vectors are three finite numbers (a size positive), an angle finite, drawing options the drawing's
  own (name, kind, range). `keepWhole` keeps atoms: with `node` or `everything` the cut is made and a warning says
  nothing was kept.
- **A twist** (`twistOf`, `twisted`): the turn about a cut's axis, from the shortest turn of its local axis onto where
  the axis points (`rotationTo`) to its rotation, measured on its local +X; which way a box's, an ellipsoid's or an
  elliptic cylinder's other axes point. `twist` sets it after `normal`. It jumps where the axis points against the local
  one (the shortest turn's axis flips there): a value shown and set, not a direction's turn.
- **A new direction** (`normal`) turns the cut from where its axis points by the shortest turn (`turnedTo`): continuous
  as the direction moves, small edits turn the other axes little; not by keeping the twist, which would spin a box about
  its axis near that point.
- **A new cut** is placed by its target's atoms (`placeCut`): a sphere or a box around a small target, with room; a
  window toward the viewer into a large one (one centered inside would show nothing); a plane through the center facing
  the camera; a cylinder a shaft along the view. With nothing selected, everything drawn.
- **Add again** returns the cut like it on those representations (shape, flip, place, direction, size within 0.5 Å and
  ~2.5°: `nearlyEqual`) instead of stacking another.
- **Another shape** keeps the place, the direction (a plane's normal is a cylinder's axis, a cone's opening) and the
  twist, the size (a box's, an ellipsoid's and a cylinder's three sizes each other's) and the flip (`reshape`); a plane
  or cone has none, so the new one is sized from what it cuts.
- **Directions** are unit; one of length 0 is refused.
- **A refused edit** reloads the form from the cut (no state change would).
- **Settings** (`settings`): `variant` (pixel: through, pixel by pixel; instance: whole copies left out by their center)
  and `solidSurface` (a cut surface closed by a cap). The variant is one for the scene: the component options' clip and
  every cut representation's (their own clip, which the options leave), and a cut written takes it (`withObject`). Why:
  Mol\* keeps it per representation, so a per-cut choice would change the other cuts on the same representations. Solid
  surface is the look's default (api.md §2, "Look"): a representation with its own keeps it. `ghost`: what the cuts
  remove shown as a faint ghost, the canvas's (`Canvas3DParams.clipGhost`, saved with its props); refused without a
  canvas. Its look (`ghostColor`, `ghostOpacity`, `ghostEdgeOpacity`, `ghostResolution`: the canvas's `clipGhost`
  params) is the canvas's while on; while off the area keeps it and turning the ghost on takes it. Why: the canvas's off
  state holds no params, and a look set and then lost on off is surprising. What neither the canvas nor the options hold
  (that look, `fadeNormal`) is the area's custom state (`prism-viewer.clipping`, `registerCustomState`): checkpoints and
  sessions carry it, so Undo of a `clipping.settings` step and a ghost turned off by Undo keep it as it was. A
  representation can opt out (its `clipGhost` param).

## 3. The section

- **Cuts · + Add**: a plane through what the selection is drawn as (narrowed by Through), or everything drawn; the new
  cut (or the one like it) becomes active. Keep the selection whole (§5) for new cuts.
- **One active cut**: the one whose editor is open (its name toggles it; Add and Clip this… make theirs), its handles
  always shown. Why: a cut drawn in the scene is never one that cannot be dragged.
- **A row**: the shape's icon and name, what it cuts (a count past three), drawn in the scene (the `scan` toggle),
  remove.
- **The editor**: the shape (a strip, as a representation's type), cut away the other side, then folds closed at first:
  Cuts through (x/N: by type, then each, those it is on first, searched past eight), Kept whole, Position and size, In
  the scene (while drawn). A box and an ellipsoid have their up, twist and three sizes; a cylinder its axis, twist, two
  diameters and height. Forms send only the fields changed, one update in flight.
- **Clip…** in the selection menu (`ui.selectionActions`) opens the section; **Clip this…** in a row's menu
  (`ui.nodeActions`: components, representations, volumes, particles) adds a plane through what the row draws, active,
  the selection left as it is.

- **Settings**: a fold under the cuts, the variant (Cut by), Solid surface, Fade arrow and Ghost, then the ghost's look
  while it is on (color, opacity, edges, resolution); its meta says "whole copies" while cuts cut by instance.

## 4. Gizmos

- **On the gizmo library** (`@molstar/graphics/util/gizmo`): its controller draws, hovers, drags, previews, commits and
  cancels; clipping is a binding (`sources/clipping/gizmo/binding.ts`: a cut as a pose, previews on the representations'
  clip, a commit through `clipping.update` by the cut's object, what a drag is on changing in the state invalidates it)
  and its layouts (`geometry.ts`). Shown for the active cut while the section is (`onActive`); their clicks claimed
  (`InteractivityManager.clickOwners`: no focus, no select).
- **Each handle says what a drag does** (the library's `Handle`): move (the screen's plane), slide (along an axis),
  plane (in a plane), rotate (about an axis), roll (about the view axis), turn (the screen's axes, as the trackball),
  scale (whole, a box across its edges, an ellipsoid across a section, a cylinder's rims about their centers: by the
  cursor's distance from the pivot over the grab's), faces (a box's face knob, an ellipsoid's tip, a cylinder's side or
  cap knob: the opposite one stays), and the cone's angle (the binding's own: its rim under the cursor). Outlines are
  handles.
- **The arrows are always drawn**, also pointing at the viewer (`clipping.settings` `fadeNormal`, off, fades them as the
  move handles' arrows do; the area's custom state, as the ghost's look). Seen end-on the center is picked over a head
  (its hit disk in front of every part).
- **Bodies** are the clip volumes' walls (`graphics/geo/clip-volumes.ts`), drawn on both sides as fronts; handles glow.
- **Above the scene** (the controller's overlay layer, out of the scene's bounds and its clip): the walls x-ray where
  the structure is in front of them, the handles on top, so a handle inside the structure can be grabbed (the library's
  defaults, `DefaultOccludedAlpha`).
- **In two parts** (`gizmo/geometry.ts`): the shape at the scene's scale, the move handles at a fixed size on screen.
  The shape is the clip object and the handles that size it (a box's edges, corners and face knobs; an ellipsoid's the
  box's on its surface: its three sections through two axes, its six tips, a knob at each octant's middle sizing it
  whole; a cylinder's rims, its four side knobs halfway up and a knob on each cap; a cone's base rim), where the clip
  is: a cut is a place in the scene; each part a set size on screen, the edges and rims as lines, the knobs as round
  points. The move and turn handles are the library's transform handles (`transformHandles`, Blender's move and rotate
  gizmo, as the move handles'), their groups first: about the cut's origin at 90 CSS pixels (`MoveScreenSize`: parts of
  a screen size of their own, `GizmoPart`), so zoomed out on a large structure they are not a few pixels, zoomed in they
  do not hide the cut. By cut what changes it (`moveOptions`): a plane's center slides along its normal, its one arrow
  points to the side it cuts away, no squares, its rings about the two axes in it; a round cylinder and a cone no ring
  about their own axis (round about it, a turn changes nothing); an elliptic cylinder (radii x and z apart), a box and
  an ellipsoid all. The move handles have hit parts (thick meshes); the shape's lines and points are a pick pixel at
  least. A turn hides the other handles, a move the turning and sizing ones (the controller's).

## 5. Keep whole

- **Mol\*'s clipping layer** on a structure representation (the cell `representation.clipping` writes, tagged
  `clipping-controls`): group bit _i_ means "not clipped by clip object _i_". A cut keeps atoms by its object's index
  there. Only the first six (Mol\*'s groups).
- **Bits on the atoms' flags**: an atom's flags are those of the last layer covering it, so a bit is set or cleared on
  them (`withKept`, `withoutKept`), never by a layer over them.
- **Removal remaps** in the same step: a cut taken off a representation drops its bit and the later bits move down
  (`afterRemoval`); the component manager replaces its own objects in place (`replaceOwned`), fewer of them leaving a
  `none` object (clipping nothing, never a cut) to hold a place before one of the representation's own. Why: an index is
  what a group means.

## 6. Drawn in the scene

- **`ClipVolumesShape3D`** (plugin core, imported by the feature: registered on import, Prism only): from a
  representation that can be cut to a shape provider of its clip objects (by tag) as translucent walls and solid edges
  (`ShapeProvider.transparency`, per group); options for walls, edges, colors, opacity, width, extent (a plane's disk
  and a cone's reach relative to what they cut, 0.25–4×: in proportion when that changes). Built from its data (objects,
  bounds, options): the shape representation keeps its first `getShape`.
- **One drawing a cut**, deduplicated by its tag (`clip-volume-<cut tag>`; every cell also `clip-volumes`), under the
  cut's first representation; it moves when a `clipping.*` run takes the cut off that one, and goes with the cut. Saved
  with the state. A representation removed otherwise takes its drawing with it (a child): the cut shows as not drawn,
  and the scan toggle draws it again. Why: redrawing it outside a step would leave a second one when the removal is
  undone.
