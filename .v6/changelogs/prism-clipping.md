## Prism: clipping

- Prism: a clipping feature in the default spec. The Selection tab's Clipping section adds cuts (plane, sphere, box,
  cylinder, cone) through what the selection is drawn as (or everything drawn), narrowed by representation type; the
  shape is picked on the cut and changed in place. Clip… in the selection menu, Clip this… in the scene tree's row
  menus (components, representations, volumes, particles) adds a cut through what the row draws.
- Prism: one cut is active at a time, its handles shown in 3D: move, slide, turn, size; a box by its face knobs and its
  edges, a cylinder by its rims, a cone by its base rim; a drag runs live and is one undo step; Esc cancels.
- Prism: a cut can leave atoms whole (`clipping.keep`, Mol*'s clipping layer), be turned on or off per representation
  (by type, one step), and be drawn in the scene (`clipping.show`, saved with it). The `clipping` area and skill:
  `list`, `add`, `update`, `apply`, `keep`, `release`, `show`, `remove`, `clear`.
- `Clip` objects have an optional `tag` (not used to clip): a clip manager's cut is the objects sharing its tag
  (`RepresentationClip.sameCut`).
- `ClipVolumesShape3D` draws a representation's clip objects as a shape (translucent walls, solid edges);
  `ShapeProvider.transparency` gives a shape per-group transparency (applied by `ShapeRepresentation3D`).
- `VisualSet` (`@molstar/graphics/util/visual-set`) holds representations drawn outside the state tree (gizmos, guides,
  shadow molecules); `PluginContext.clickOwners` lets a tool own its representations' clicks.
- The component manager's clip objects (its options) are replaced in place, the representations' own kept; fewer of
  them leave a `none` object holding a place, so a representation's own objects keep their indices.
- A representation moved by its transform keeps its bounds (the camera's near and far planes, culling, picking).
- Prism: `ui.nodeActions` adds rows to the scene tree's node menus; the conformance checks hold for extensions' Prism
  code too (no rxjs in views, scoped attachments, owned subscriptions, the plugin only through operations, a command
  for every operation); the canvas' bottom bar no longer blocks picking beside its items; a trackpad pinch never zooms
  the page.
- Prism: the `clipping` operations take a cut by its tag too (the views and gizmos use it); vectors, angles and drawing
  options are checked; a box keeps its twist when its direction changes; a refused edit reloads the form.
