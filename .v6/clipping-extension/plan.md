# Clipping extension: plan (working notes, deleted before the PR)

Before the PR, what holds moves to `extensions/clipping/_spec/`, and this folder is deleted.

**Prism particles are done (2026-10-21, viewer-spec.md §3.2)**: they are an opt-in registry entry with `particles` and `shape` node kinds. The selection menu is done (viewer-spec.md §3: `ui.selectionActions`, where this extension's Clip… goes; the tree's selection filter). This plan is next. Steps 1 and 2 here are done.

## Goal

`extensions/clipping` (`@molstar/clipping-extension`) brings:

1. **A Clipping section in Prism's Selection tab.** It is a sibling of "Color, transparency, glow". It lists the cuts, and adds, edits and removes them. It applies to every **clippable representation**:
   - structures, volumes, shapes (MVS primitives, meshes, PLY) and particles;
   - that is, any representation cell whose geometry carries `clip`.
   - The current selection finds the representations a cut applies to. A filter by representation type narrows that further: a selection drawn as ball-and-stick and as a surface can be cut on the surface only.
2. **Gizmos** in the 3D view for each cut, of any shape the data model has: plane, sphere, cube, cylinder, infinite cone. They are toggled from the section. They are raw visuals added to Canvas3D, not state cells. Dragging them moves, turns and sizes the cut live, and Esc cancels.
3. **"Clip this…"** in the row menu of any clippable node in the scene tree. It opens the Selection tab with the Clipping section unfolded, the node selected, and the type filter set to that representation. This needs a new Prism extension point, added in this PR.
4. **Exclusion** ("keep these atoms whole") in the same section: per-loci clipping layers, for structures only.
5. **`VisualSet`** in `packages/graphics/src/util/visual-set.ts`. It owns a set of raw representations in one Canvas3D: create, add, set a transform, check whether a loci is its own, and dispose all of them in one call.
6. **Storybook stories** for the scenarios below.

## Decisions (David, 2026-10-21)

- **D1. The state tree is the source of truth.**
  - The cuts are the clip objects in the representation cells' params. The extension keeps no store of its own.
  - A clip set from code (`type.params.clip` through a state update, an MVS `clip` node, mesoscale-explorer's group clip) shows in the section too.
  - So the section's list is *derived* from the state. Gizmos are a view of that list.
- **D2. Clipping is per representation.**
  - The selection only finds representations; it does not narrow the cut to atoms.
  - It has to work for every object whose representation has `clip`: structures, volumes, shapes and particles.
- **D3. Live drag updates the live representations outside the state, then makes one state update on release.**
  - Be careful with cancels. See *Drag and cancel*.
- **D4. The scene tree row-menu point is added in this PR.**
- **D5. "Exclude from clipping"** (`representation.clipping`, the layers by loci) is shown in the section too.
- **D6. Every shape the data model has:** plane, sphere, cube, cylinder, infiniteCone. Each shape gets a gizmo.

## What exists already

**Clip objects:**
- `Clip.Props` (`core/util/clip.ts`) is `{ variant: 'instance' | 'pixel', objects: { type, invert, position, rotation: {axis, angle°}, scale, transform }[] }`, with at most 16 objects.
- It is part of `BaseGeometry.Params` (`graphics/geo/geometry/base.ts:154`), so every geometry has it.
- It sits in a cell's params at:
  - `type.params.clip` for `StructureRepresentation3D` and `VolumeRepresentation3D`, and probably particles (to check);
  - `clip` at the top for `ShapeRepresentation3D`, whose params are the provider's own geometry params.
- `clipPrimitive` (spheres only) clips whole spheres.
- It is tested in world space, after the representation's transform. `variant: 'instance'` tests the instance center, which mesoscale uses for big instanced sets.
- It is not part of `Representation.State`, so a change goes through `createOrUpdate`. A clip-only change should touch only uniforms; to be measured.
- Each shape's local frame is documented in `extensions/debug-helpers/src/clip-object-helper.ts`:
  - plane normal: rotated +Y;
  - sphere and cube: `scale * 0.5`;
  - cylinder: axis Y;
  - cone: axis Z.

**Exclusion layers** (the `Clipping` theme, `setStructureClipping`, cell tag `clipping-controls`):
- Structures only, by loci.
- Bit *i* means "not clipped by object *i*", for objects 0–5 only.
- Prism has `representation.clipping`, but no UI for it. Its layers show in the applied-layers list and can only be removed there.
- There is no `clearStructureClipping`.

**The component manager clobbers clip.** It writes `options.clipObjects` into every structure representation's `clip` whenever its options are applied (`plugin/core/state/manager/structure/component.ts:152-197`). Prism exposes this as `component.options({ clip })`. Under D1 this would wipe the cuts, so the manager has to change; see *Plugin core changes*.

**mesoscale-explorer** (`apps/mesoscale-explorer/src/data/state.ts:96-205`, `ui/entities.tsx:975-1000, 1499-1515`) is the one app that clips today:
- One clip object per group or entity, edited in `SimpleClipParams`, relative to the scene's bounding sphere. Position and scale are given as % of the radius.
- It is written to every filtered entity's representation, in `type.params.clip.objects` or `clip.objects` for shapes, which are the two param paths. It is mirrored into the group cell's own `clip` param, so the UI can read it back.
- It picks `variant: 'instance'` for big instanced sets and `'pixel'` for single ones, and sets `clipPrimitive: true` for spheres.
- Lessons:
  - The two param paths need one adapter.
  - Relative units confuse when the scene changes. We store absolute Å in the state, and show relative helpers only in the UI if at all.
  - The variant is a per-representation choice, so a cut must not override it.
  - It replaces `objects` wholesale, so a second cut from our section would be lost by its next edit. We should co-exist: our writes keep objects we do not own.

**Selection to representations:** `resolveCovered` / `applyLayer` (`prism/viewer/src/api/representation.ts:780-925`) work for structures. Volumes and shapes need their own resolution: the selection's structures, plus picked or selected scene nodes. See *Targets*.

**Prism extension points:**
- `ui.sections` with `slot: 'selection'` (Measure uses it), plus `api` areas, commands and skills. The area name `clipping` is free.
- Row menus are built by the viewer (`SceneNode.menu`, `views/structure/scene-tree.tsx:1147`). We add the point.

**Canvas3D public API**, enough for gizmos:
- `add` / `update` / `remove`, `getLoci` through `repr.getLoci`, `mark`, `interaction.hover`, `input.drag` / `interactionEnd`.
- Trackball bindings can be swapped through `setAttribs({ trackball: { bindings } })`.
- Raw representations survive a WebGL/WebGPU switch.
- No Canvas3D or gl changes (WebGPU freeze). `handle-helper.ts` is the model for group ids and highlight-only marks.

## Lessons from stories-untold `pmg_anim`

From `designs/clip-guides`:
- **Pickable parts must be opaque.** Picking needs alpha ≥ 0.5, so the translucent faces are their own unpickable visual, and the handles are opaque.
- **Draw a plane's disk where the plane crosses its targets**, not at `position`, sized to their reach. The same goes for a cone's and a cylinder's extent: they are infinite or long, and are drawn only over the targets.
- **Live drag cannot go through the state.** It lags. They wrote `uClipObject*` uniforms by hand; we do D3 instead.
- **MVS plane to rotation is singular for a normal pointing down (H005).** It is fixed in its own session (`mvs-clip-plane-with-down`). Our plane gizmo's normal-to-rotation conversion must handle the same case.
- **Gestures:**
  - Move: along the screen axes, Å per pixel.
  - Rotate: 0.5° per px.
  - Size: ×2 every 150 px.
  - Flip: toggles `invert`.
- **New cuts are placed** at the targets' bounds:
  - A plane goes through the center, facing the camera.
  - Closed shapes go at 0.85 r on the near side. One centered in a hollow capsid cuts nothing.

## Design

### A cut is derived, not stored (D1)

**Clip adapters.** A small registry maps a representation transformer to `{ get(params): Clip.Props, set(params, clip) }`:
- built in for structure, volume, shape and particles;
- open for extensions (`ClipAdapters.register`);
- **Where it lives:** `packages/plugin/core`, so the component manager and the section share it. Or in the extension, if the manager fix doesn't need it.

**Reading the state.** The extension scans the clippable cells (the adapter knows them) and groups their clip objects into **cuts**:
- Objects that are equal (`type`, `invert`, `position`, `rotation`, `scale`, `transform`; `Clip` already has equality helpers) on several cells make one cut, applied to those cells.
- `cut = { key, object: Clip.Object, targets: StateTransform.Ref[] }`, where the key is the object's signature.
- Two cuts with the same geometry on different targets merge. That is acceptable: it is the same cut.
- This needs no id in the params. An optional hidden `label` on a clip object, for naming cuts, would be a core change; deferred.

**Editing a cut** rewrites the matching object on each of its target cells, in one state update, through the adapter. Everything else in that cell's `clip` is kept:
- objects we do not own (mesoscale's, the manager's global ones);
- the `variant`.

**Order matters for exclusion layers** (bit *i* is object *i*). Adding a cut appends it; removing one shifts the indices of the rest. The section therefore re-maps the exclusion layers of the cells it touches, or warns that removal shifts them. It is decided in step 5 with a test.

### Targets

`clipping.add({ target, types? })` resolves:
- **Structures:** the target's components → representations covering it (the existing `resolveCovered`), then filtered by type name.
- **Any clippable cell:** `node?: NodeParam()` (Prism's node selectors, `_spec/nodes.md`: `volume > representation`, `#ref`, `structure[label="1CBS"] component[type=ligand] > representation`). This covers volumes and shapes. Representation nodes under a matched node are taken, and other kinds are skipped. "Clip this…" passes `#<ref>`.
- **Both given:** the union. **Neither given:** `target: 'selection'`.
- **`types` options:** the representation types drawn in the resolved set, for example `cartoon`, `molecular-surface`, `isosurface`, `mesh`.

### API area `clipping` (headless, `prism-viewer.ts`)

Operations, each also a command and in the skill:
- `clipping.list`: the derived cuts, each with its targets and the exclusions in its structures.
- `clipping.add`: `{ shape, target = 'selection', types?, invert? }`, placed from the targets' bounds and the camera.
- `clipping.update`: `{ cut, position?, rotation?, scale?, invert? }`, plus a target change (add or remove representations).
- `clipping.remove`: `{ cut }`. `clipping.clear`: `{ target? }`.
- `clipping.exclude`: wraps or replaces `representation.clipping` for D5, with the cut picked by index, at most 6.
- `clipping.gizmos`: `{ show, cut? }`.

Everything runs as the person, in one `dataTransaction`. It returns `already` when nothing changed.

### UI (`prism-viewer-ui.tsx`, lazy)

A section with `id: 'clipping.section.clip'`, `slot: 'selection'`, `icon: 'crop'`, `order: 0.5`. It contains:
- **Add:** a shape picker over all five shapes, plus the type filter, a multi-select of the types drawn in the selection, "All" by default.
- **The list of cuts.** Each row has:
  - a shape icon, the target count (say "3 representations"), flip, a gizmo toggle and remove;
  - unfolded: numeric fields (center, scale, normal or axis and angle) and the targets, which can be removed one by one.
- **Keep whole:** the exclusion for the selection, per cut (first 6).
- **"Show gizmos"** for all cuts.

Prism rules apply: the DOM runtime, no rxjs in views, `scoped` attachments, `--msp-*` tokens, and checks in light, dark, wide and narrow.

### Gizmos (`gizmos.ts`, on `VisualSet`)

**Per cut, two `ShapeRepresentation`s.** Both get `markerActions: Highlighting`, `clip: none` and `pickable` set as listed.

**Faces** (translucent, unpickable), by shape:

| Shape | Faces |
| --- | --- |
| Plane | a disk over the targets' reach |
| Sphere | a UV sphere |
| Cube | 6 faces |
| Cylinder | a capped tube |
| Infinite cone | a cone over the targets' reach |

**Handles** (opaque, pickable). The groups are the parts, and each label names its action:

| Part | What it does |
| --- | --- |
| Center knob | moves the cut in the screen plane |
| Normal or axis arrow | slides along it |
| Rings | turn the cut (plane, cube, cylinder, cone) |
| Face or end knobs | size the cut along an axis |
| A flip button-knob, or the arrow's head | toggles `invert` (Flip is in the section as well) |

**Built once in unit space**, moved by `setState({ transform })` and `canvas3d.update`. A rebuild happens only when the aspect changes.

**Side effects:**
- **A click on a gizmo must not focus the camera or change the selection.** It needs a plugin opt-out; see below. Check Prism's own click path as well.
- **Camera reset:** keep the gizmos inside the targets' sphere, or `manualReset` while they are shown.
- **Gizmos are not in the session or the snapshot.** The gizmo set follows the derived cuts: on each state change it adds, updates or removes gizmos.

### Drag and cancel (D3)

**On press, over a handle group:**
- Snapshot the full trackball bindings, then empty the drag bindings.
- Record each target cell's ref, its current `clip` value from the *state params*, and its state version.

**On move:**
- Compute the new object.
- Call `repr.createOrUpdate({ ...repr.props, clip })` on each live representation (`cell.obj.data.repr`), outside the state.
- Move the gizmo transform.
- Coalesce the work to one update per frame.

**On release:** one state update to the new object, with an undo step. Restore the bindings.

**A cancel restores the live representations from the recorded state `clip`.** Bindings are restored in every case, in a `finally`. A cancel comes from:
- Esc, while the canvas or document has focus;
- a lost release: pointer cancel, blur, `interactionEnd` without a drag end;
- the target cell changing or being removed during the drag (a state event on one of the recorded refs: undo, a preset, another edit). That is a cancel, not a commit over it;
- the gizmos being turned off, or the extension unmounting.

**The release commits only if the cells are still at the recorded versions.** Otherwise it cancels.

**Tests:**
- The pure gesture math.
- The cancel paths, against a headless plugin with a fake input:
  - Esc;
  - state changed mid-drag;
  - removed mid-drag;
  - the release lost.
- Each test asserts the live props equal the state params afterward.

### `VisualSet` (`packages/graphics/src/util/visual-set.ts`; done in step 2)

**What it holds:** any representations drawn outside the state tree: gizmos, guides, and *shadow molecules* an editor previews. It handles structure, volume and shape representations alike.

**Its API:**
- `add(repr)`: an already built representation.
- `create(build)`: a representation built with the canvas's graphics.
- `addShape(getShape, utils, data)`.
- `fromProvider(provider, data, { themes, themeData, properties })`: themed as the state tree would theme it. Custom properties and theme dependencies are attached, and released when it is removed.
- `update`, `setState`, `setTransform`.
- `has(repr)`, `owns(eventLoci)`: an interaction event's `repr`, so a shadow molecule sharing a `Structure` with the tree is told apart.
- `remove`, `clear`, `dispose`.

**How it behaves:**
- Auto-synced (`syncManually: false`).
- A build that finishes after `clear` is destroyed and not added.
- It is tested with a fake canvas (`util/_test/visual-set.test.ts`).

**Not handled here:** the plugin behaviors acting on clicks (focus, selection) on its representations. That is the opt-out in plugin core. It applies to shadow molecules too, whose `StructureElement` loci must not reach the selection manager.

### "Clip this…" and the row-menu point (D4)

**The point**, in `packages/prism/viewer`:
- A headless `ViewerExtension.nodeMenu?: { id, when(node: SceneNode): boolean, item(node, ctx): SceneMenuItem }[]`.
- The tree appends these items to `menuOf(node)`, before Remove.
- It is recorded in `viewer-spec.md §3`, and the stub types (`LazyViewerExtensionStub.parts`) are checked as well.

**The item** shows for any node whose representations have a clip adapter: a component, a representation, a volume or a shape. It:
1. selects the node;
2. sets the section's type filter to that representation, for a representation row;
3. opens `view.open({ tab: 'selection', section: 'clipping.section.clip' })`.

### Plugin core changes

1. **The component manager stops clobbering `clip`.**
   - `clipObjects` becomes the *global* objects. On an options change it replaces only the objects equal to the previous global ones, appended at the end, and keeps the rest.
   - This is tested.
   - It goes in `breaking-v6-changes.md` if the meaning changes for callers.
   - The focus representation behavior keeps getting only the global ones.
2. **The clip adapter registry**, if it lives in plugin core.
3. **The gizmo click opt-out.** `FocusLoci` and `SelectLoci` skip loci from a representation marked as a helper: for example a repr state flag, or a `ShapeGroup` shape whose source data is tagged. The mechanism is decided in step 3.

## Steps

**Step 1 is done (2026-10-21).** A clip-only `repr.createOrUpdate({ ...props, clip })`:
- **makes no render object anew and touches no attribute, texture or non-clip define**, on 1cbs cartoon, ball-and-stick, molecular surface and gaussian surface, a test volume isosurface and a shape mesh;
- **takes 0.01–0.12 ms an update**, outside GL;
- **changes only uniforms.** Besides the `uClipObject*` uniforms and `dClipObjectCount`, it re-sets `uLod`, `uInteriorColor` and `uInteriorSubstance` every update, though their values stay the same. That is cheap.
- **The guard** is `packages/graphics/src/repr/_test/clip-update.test.ts`.
- **Consequence:** the live drag updates each target's live representation every frame; no uniform hacks.
- **Not measured yet:** the release's state update. It also re-creates the theme (`StructureRepresentation3D.update` calls `setTheme`), to check in step 5 with a headless plugin. It happens once per drag.

1. **Spike, measuring D3.** A clip-only `createOrUpdate` on a live representation, for a structure cartoon, a molecular surface on 1cbs, a volume isosurface and a shape mesh: does anything rebuild, and how long does it take?
2. **`VisualSet` and its tests.**
3. **Plugin core (done 2026-10-21):**
   - `RepresentationClip` (`state/helpers/representation-clip.ts`) reads and writes a cell's clip at either path, so no adapter registry is needed: every representation transform holds it at `type.params.clip` or `clip`, and extensions' do too.
   - The component manager replaces only its own objects (`replaceOwned`).
   - `PluginContext.clickOwners`: a tool claims its representations' clicks, which then never reach `behaviors.interaction.click` (focus, select, host events). Hover is untouched.
4. **`extensions/clipping` package wiring:** `inventory.json`, root tsconfig refs, `apps/prism-viewer` deps and `extensions.ts`.
5. **The `clipping` API area:** derive the cuts, the operations, commands and a skill, and exclusion re-indexing. Quick-lane tests through a headless viewer, including a clip set from code showing up in `list`.
6. **The Selection section UI.**
7. **Gizmos for the five shapes:** hover, drag, commit and cancel, with the gesture and cancel tests.
8. **The row-menu point in Prism and "Clip this…".**
9. **Stories** in `apps/prism-viewer/src/_stories/clipping.stories.tsx`, with fixtures from `ExampleFiles` and no network:
   1. a plane on 1cbs, surface only, while ball-and-stick stays whole;
   2. a sphere window;
   3. a cube;
   4. a cylinder and a cone;
   5. a volume isosurface cut (local fixture to find);
   6. a shape or MVS primitive cut;
   7. a clip set from code, listed in the section;
   8. gizmos on, dragging a plane, then Esc;
   9. "Clip this…" from a row;
   10. "Keep whole" exclusion.
10. **Write-up:** the extension's `_spec/README.md`, a `.v6/changelogs/` entry, a `CHANGELOG.md` line, then this folder is deleted.

## Open questions

- **Answered:** volumes and shapes are reached through node selectors (`node`), not the atom `target`. Prism's selection holds atoms only.
- Should a node that is picked in the tree (the inspector's subject) count as "the selection" for the section when no atoms are selected? Proposal: yes, the section shows "Applies to: <node>".
- Particles: where does its `clip` live (`type.params.clip`?), and does Prism show particles at all?
- Should a cut be able to change a representation's `variant`, as mesoscale picks `instance` for big instanced sets? Proposal: no. Keep the cell's variant, and show it read-only.
