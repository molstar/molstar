## Trajectories: molecules whole across the periodic box

- Frames of a simulation in a periodic box show every molecule whole: a molecule the engine wrapped across a face is
  no longer drawn with bonds and backbone stretched across the box. On by default for frames with a simulation box
  (coordinates with a cell: XTC, TRR, DCD, NetCDF; GRO); `ModelFromTrajectory`'s new `periodic: 'off'` shows the frames
  as written. Crystal and NMR models are never changed.
- `periodic: 'box'` also moves each molecule into the drawn box, and `periodicCenter` centers the box on a selection
  (polymers by default): the protein in the middle, the solvent around it, in every frame.
- GRO reads triclinic boxes (all nine box values: the rhombic dodecahedron, the truncated octahedron); it read a
  rectangular box before. `GroHeader.boxVectors` holds them.
- Library: `PeriodicBox` (`@molstar/core/math/geometry/periodic-box`), `Model.PeriodicBox`,
  `Model.withAtomicConformation`, `ModelPeriodicity` (`makeWhole`, `apply`) (`@molstar/model/model/structure/model/periodic`).
- Prism stories: `Viewer › Simulations` (Periodic box, Trajectories per format, Simularium moved there).
