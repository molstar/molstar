# Moving models by hand

The `move` feature (`src/features/move`; opt-in: not in the default spec, the Prism app adds it): the Selection tab's
Move section and the gizmo library's transform handles (`@molstar/graphics/util/gizmo/transform`) on a model: a
structure, a particle list, a shape or a volume.

- **Through the model's transform decorator**: a structure's conformation transform (`structure.transform`, Mol\*'s
  `TransformStructureConformation`), a particle list's (`particles.transform`, `TransformParticleList`: the coordinates
  moved, the orientations turned, given where it turns a list without them), a shape's (`shape.transform`,
  `TransformShapeProvider`: its instances' transforms premultiplied), a volume's (`volume.transform`, Mol\*'s
  `VolumeTransform`: its grid placed); no state of its own. A rigid move is a turn about the model's center and a move,
  `T(position) R T(−center)`: the handles' frame is the center as placed now, turned as now (`placementPose`,
  `poseMatrix`). A decorator is inserted on the model as the hierarchy has it (its data the model as loaded, the
  representations under the decorator); `sources/particle-transform.ts`, `sources/shape-transform.ts` and
  `sources/volume-transform.ts` read them.
- **Each kind's center**: a structure's boundary, a list's, a volume's grid; a shape's is known only once it is made:
  the bounds of its representations as drawn, all its instances (`shapeDrawnCenter`; handles once it is drawn).
  `shape.transform` turns about the same center (`shapePivot`), so a turn by the API and one by hand agree; not drawn,
  about the shape made once with its default props, all its instances (`getShapeProviderCenter`). It writes a matrix,
  which the handles read at once.
- **What can be moved**: the structures, the particle lists, the shapes loaded on their own (`sources/shapes.ts`: not a
  structure's, a map's or particles'; a shape's representations found under its transform), the volumes loaded on their
  own (not one streamed for a structure, `isStreamedForStructure`: it moves with the structure; `volume.transform`
  refuses it, naming the structure's move). Why: streamed density is fetched for the structure in its frame, a region at
  a time; moved apart it no longer fits, and what is fetched next is placed by the structure, not by the move.
- **The handles at a fixed size on screen** (the gizmo layout's `screenSize`, in the overlay layer, on top of the
  model): a model's size and the zoom don't change them, so a large one's are not off screen and a small one's are not
  too small to grab.
- **Blender's move and rotate gizmo** (`transformHandles`' defaults, no scale): the center moves in the screen's plane,
  an arrow head slides along its axis, a small square moves in its plane; a ring (its front half) turns about its axis
  by the angle the cursor goes around the center on the screen (Blender's; `RingRule`, in `gestures.ts`), its radius at
  the start and now drawn; the grey circle around them rolls it about the view axis; the area inside the rings turns it
  freely, as the trackball turns the scene. It takes the camera's drags over that area: moving a model is a mode entered
  on purpose (the section shown).
- **A drag's preview moves the model's representations as drawn** (their transient transform, `previewTransform`:
  applied on the GPU, no instance, grid or culling layout touched); recomputing it each frame would not keep up, nor
  would baking a transform into a mesoscale model's instances. Only the representations of its kind: a volume streamed
  under a structure is not moved by the transform, nor a list's unit cell (centered on the particles again, unturned: a
  periodic box is the simulation's axes). Its release is one transform operation with the matrix from where it was
  loaded (a step); the transient transform is put back once the state has it, the canvas's frames held in between
  (`Canvas3D.holdFrames`, a few seconds at most: the representations rebuilt still carrying the preview would show the
  model moved twice). A change of the model's cells under a drag cancels it (a list's frame too), as does the model
  going away. A placement with more than a rigid move (a scale or shear given as a matrix) keeps it under the new pose
  (`applyToPlacement`).
- **A release's cost by kind**: a particle list moves what is built (its targets' instances and grids); a structure is
  built again (few instances); a shape's representations are made again (`ShapeRepresentation3D` calls the getter it was
  made with: a moved shape's is another); a volume's representations again from the placed grid (an isosurface computed
  again).
- **The section**: which model (the structure holding the selection, or the first, picked again as they come and go; a
  row each when several: the structures, the particle lists, the shapes, the volumes), Reset (the kind's
  `transform-reset`); the handles show while the section does. The `move` area adds no operations; the four transforms
  take the same arguments (`api/rigid-transform.ts`).
