# Models

`Model` (`src/model/structure/model/model.ts`): one model of a file, read-only. Atoms are the integers `[0, atomCount)`
(`ElementIndex`).

## 1. Fields

| Field                                         | What                                                                                               |
| --------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| `id`, `entryId`, `label`, `entry`, `modelNum` | identity; `modelNum` is the file's model number, or the frame index for frames from coordinates    |
| `sourceData: ModelFormat`                     | the parsed file (`kind` names the format, formats.md §1)                                           |
| `atomicHierarchy`, `atomicConformation`       | who and where: §2                                                                                  |
| `atomicRanges`                                | `polymerRanges`, `gapRanges` (atom ranges), `cyclicPolymerMap`                                     |
| `atomicChainOperatorMappinng`                 | chain → operator for models written with operators applied (sic, the name is in use)               |
| `coarseHierarchy`, `coarseConformation`       | spheres and gaussians of integrative models                                                        |
| `entities`, `sequence`, `properties`          | entity table and subtypes; sequences per entity; missing residues, component maps, `structAsymMap` |
| `customProperties`                            | descriptors, reference counts and assets of attached properties (properties.md §1)                 |
| `_staticPropertyData`, `_dynamicPropertyData` | property bags: §3                                                                                  |
| `parent`                                      | `Model.getRoot(m) = m.parent ?? m`; nothing in the repository sets it                              |

## 2. Hierarchy and conformation

- **Hierarchy** (`properties/atomic/hierarchy.ts`): `atoms`, `residues`, `chains` tables; `residueAtomSegments` and
  `chainAtomSegments` (offsets plus a per-atom index); `index` (lookups by label/auth ids, -1 when absent); `derived`
  (atomic numbers; per residue the trace, direction atoms, molecule and polymer type).
- **Sorting**: atoms are sorted by entity, asym and seq id within chains that have `label_seq_id`, so residues and
  chains are contiguous. `atomSourceIndex` (absolute rows of the file's `atom_site`) and `residueSourceIndex` map back;
  export and index-pair bonds read through them.
- **Segments**: a chain starts when entity or `label_asym_id` changes; a residue when seq id, auth seq id or insertion
  code changes, not the component, so microheterogeneity stays one residue.
- **Conformation** (`properties/atomic/conformation.ts`): `id`, `atomId`, `occupancy`, `B_iso_or_equiv`, `x/y/z`.
  `atomId` lives here because mmCIF gives corresponding atoms different ids in different models. Read positions through
  a unit's `conformation` (it applies the operator), not from `x/y/z`.
- **Coarse**: `spheres` and `gaussians` (`CoarseElements`: entity, asym, seq ranges, chain segments, polymer and gap
  ranges) with their own conformation and a lazy `CoarseIndex`.
- **Ranges** (`utils/atomic-ranges.ts`): a gap where seq ids jump or backbone atoms are more than 3 Å apart (10 Å
  coarse); computed from the first frame and shared by later ones.

## 3. Frames and what they share

A frame is `{ ...model, id, modelNum, atomicConformation, _dynamicPropertyData: fresh }`: from coordinates
(`Model.trajectoryFromModelAndCoordinates`, `trajectoryFromTopologyAndCoordinates`) and from the next model of a file
whose atom, residue and chain tables equal the previous one's (`createModels`).

| Shared by frames                                                                                 | Per frame                                              |
| ------------------------------------------------------------------------------------------------ | ------------------------------------------------------ |
| hierarchies, ranges, entities, sequence, `properties`, `customProperties`, `_staticPropertyData` | `id`, `modelNum`, conformation, `_dynamicPropertyData` |

- **Static data** does not depend on coordinates; **dynamic data** does. Neither bag is read directly: providers and
  `Model.*` helpers own their keys (`__AtomicRadii__`, `__Center__`, `TrajectoryInfo`, `CoordinatesHistory`, ...).
- `Model.areHierarchiesEqual` is identity of the hierarchies: shared topology.
- A frame must have the model's atom count (throws otherwise); frames in file order are reordered by `atomSourceIndex`,
  in place on first use (`Coordinates.getAtomicConformation`).
- `Model.CoordinatesHistory` (static, on the root model) answers whether two models have equal coordinates for a set of
  elements; `Unit.isSameConformation` asks it before comparing coordinates.

## 4. Trajectories, coordinates, topologies

- `Trajectory`: `frameCount`, `duration`, `representative`, `getFrameAtIndex(i)` (a model or a task). `ArrayTrajectory`
  holds models.
- `Coordinates` (`coordinates/coordinates.ts`): frames of `x/y/z` with optional cell, velocities, forces and a time;
  `step` when the file says which simulation step a frame is (DCD), to match box files by: `coordinatesWithXsc`
  (`formats/structure/namd-xsc.ts`) gives frames a NAMD `.xsc`/`.xst` file's boxes, the last entry at or before their
  step, the corner the file's center less half of a + b + c.
- `Topology` (`topology/topology.ts`): `BasicData` plus bonds as index pairs; combined with coordinates into a
  trajectory whose frames share one `IndexPairBonds`.
- **Periodic frames** (`model/periodic.ts`): `Model.PeriodicBox` (dynamic) is a frame's simulation box, set only by
  readers of simulation output (frames from coordinates with a cell or a `Frame.box`, GRO, LAMMPS); a crystal's cell
  (`ModelSymmetry`) is never one, so no crystal or NMR model is moved. `ModelPeriodicity.makeWhole(frame)` moves each
  atom by a lattice vector next to its parent in a plan built once per topology (a BFS forest over `IndexPairBonds`,
  else over the hierarchy's hops, the atom before in a residue and the previous polymer residue in a chain, joined
  across the molecules it leaves apart by distance links: element pairing thresholds under the minimum image, found with
  a cell list that wraps at the box's faces, never to water, ions or metals; `planFromHierarchyAndDistances`). The whole
  frame is `Model.withAtomicConformation`: new coordinates, the ranges found again on whole coordinates (shared by the
  topology's whole frames), the reader's dynamic data without `__Center__`; one model per recent source frame and mode
  (an LRU of 16 frames per topology, held weakly by its hierarchy: caching every frame would hold a trajectory twice,
  and a module-level LRU would hold a dropped one). `apply(frame, { kind: 'box', center? })` then moves each whole
  molecule by a lattice vector so its center is in the box (at the box's origin); with `center` atoms, their molecules
  are kept together, everything is translated so their mean is at the box's center, and every other molecule is moved
  into the box around them. `applyNoJump(trajectory, index)` keeps each molecule's center from moving a box length
  between frames: shifts in box vectors walked from the last frame asked or a checkpoint every 64 frames (sequential by
  nature). A frame whose file gives unwrapped coordinates beside the wrapped ones (`Frame.unwrapped`: LAMMPS image
  flags) has `Model.UnwrappedConformation` (made when first asked): whole moves each molecule's unwrapped coordinates so
  its root is where the frame has it, and no jump is frame 0's shifts on them, no walk; the frame itself stays as read.
  `Unit.remapModel` drops polymer and gap elements when the ranges differ.

## 5. Other models

- **Volume** (`model/volume/`): a `Grid` (`Tensor` cells, a spacegroup or matrix transform, stats, periodicity) and
  instances (transforms). Iso values absolute or relative (σ). Property bags `_propertyData` and `_localPropertyData`
  (caches on `volume.parent ?? volume`).
- **ParticleList** (`model/particles/`): packed coordinates, optional rotations, radii, fibers and attributes; targets
  (a structure, shape or volume per particle); lazily derived transforms and target groups.
- **Sequence** (`model/sequence/`): protein, DNA, RNA or generic sequences with seq ids and microheterogeneity;
  alignment (`alignment/`).
- **Shape** (`model/shape/`): structural only; graphics builds the renderable one.
- **Location and Loci unions** (`location.ts`, `loci.ts`): every kind the plugin can pick, label and highlight.
