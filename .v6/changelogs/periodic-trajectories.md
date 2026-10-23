## Trajectories: molecules whole across the periodic box

- Frames of a simulation in a periodic box show every molecule whole: a molecule the engine wrapped across a face is
  no longer drawn with bonds and backbone stretched across the box. On by default for frames with a simulation box
  (coordinates with a cell: XTC, TRR, DCD, NetCDF; GRO); `ModelFromTrajectory`'s new `periodic: 'off'` shows the frames
  as written. Crystal and NMR models are never changed.
- `periodic: 'box'` also moves each molecule into the drawn box, and `periodicCenter` centers the box on a selection
  (polymers by default): the protein in the middle, the solvent around it, in every frame.
- `periodic: 'nojump'` keeps each molecule from jumping a box length between frames (molecules leave the box as they
  diffuse): seeks walk from a checkpoint every 64 frames, playing steps one frame.
- Prism: `trajectory.periodic` (whole, in the box, centered on the polymers, no jumps, as in the file) and its command
  ("periodic box", "pbc", "unwrap", "center protein in box"); the assistant's `periodic_box` tool; a box menu at the end
  of the trajectory bar for a structure in a simulation box.
- MolViewSpec: `custom.molstar_periodic` on a structure node (`'whole' | 'box' | 'nojump' | 'off'`, or
  `{ mode, center }` with center `'box'`, `'polymers'` or a selection expression), a Mol\*-specific key.
- The ensemble structure of all frames (`StructureFromTrajectory`) takes the same `periodic` params: whole by default.
- LAMMPS: triclinic boxes (tilt factors, the dump's bounding-box bounds corrected), the box's corner as its origin
  (the unit cell is drawn there), scaled coordinates through the tilt; frames get a box only when periodic in x, y, z.
- GRO reads triclinic boxes (all nine box values: the rhombic dodecahedron, the truncated octahedron); it read a
  rectangular box before. `GroHeader.boxVectors` holds them.
- Library: `PeriodicBox` (`@molstar/core/math/geometry/periodic-box`), `Model.PeriodicBox`,
  `Model.withAtomicConformation`, `Frame.box`, `ModelPeriodicity` (`makeWhole`, `apply`, `applyNoJump`) (`@molstar/model/model/structure/model/periodic`).
- Prism stories: `Viewer › Simulations` (Periodic box, Trajectories per format, Simularium moved there).
