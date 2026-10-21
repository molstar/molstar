# Moving structures by hand

The `structure-move` feature (`src/features/structure-move`; opt-in: not in the default spec, the Prism app adds it):
the Selection tab's Move section and the gizmo library's transform handles (`@molstar/graphics/util/gizmo/transform`) on
a structure.

- **Through the conformation transform** (`structure.transform`, Mol\*'s `TransformStructureConformation`): no state of
  its own. A rigid move is a turn about the structure's center and a move, `T(position) R T(−center)`: the handles'
  frame is the center as placed now, turned as now (`placementPose`, `poseMatrix`).
- **A drag's preview moves the structure's representations as drawn** (their state transforms, on top of what they had);
  recomputing a structure each frame would not keep up. Only structure representations: a volume streamed under the
  structure is not moved by the transform. Its release is one `structure.transform` with the matrix from where it was
  loaded (a step); the representations' own transforms are put back once the state has it, the canvas's frames held in
  between (`Canvas3D.holdFrames`, a few seconds at most: the representations rebuilt still carrying the preview would
  show the structure moved twice). A change of the structure's cells under a drag cancels it, as does the structure
  going away. A placement with more than a rigid move (a scale or shear given as a matrix) keeps it under the new pose
  (`applyToPlacement`).
- **The section**: which structure (the one holding the selection, or the first, picked again as structures come and go;
  a row each when several), Reset (`structure.transform-reset`); the handles show while the section does. The `move`
  area adds no operations.
