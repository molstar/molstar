# Plan

Goal: a frame from a periodic run shows every molecule whole, and optionally wrapped around a chosen center, for every
trajectory source (coordinates + topology, multi-model files, LAMMPS), in the library, the plugin, Prism and headless.
Non-goals for now: bonds and interactions between periodic images, drawing neighbouring images (a later P1 supercell
option, §7), reading NAMD `.xsc` files.

## 1. Operations

Three, as in every MD tool ([analysis §5](analysis.md#5-how-other-tools-do-it)), each a pure function of a frame's
coordinates:

| Operation   | Does                                                                                     | Needs                                       | Per frame |
| ----------- | ---------------------------------------------------------------------------------------- | ------------------------------------------- | --------- |
| **whole**   | each atom moved by a lattice vector to the minimum image of its parent in a spanning tree | the lattice (no origin), a tree, once       | O(N)      |
| **wrap**    | each molecule moved by a lattice vector so its center lies within half a box of a center | whole molecules, a center (box or selection) | O(N)      |
| **no jump** | each molecule's lattice shift chosen so its center moves less than half a box per frame  | the previous frame's shifts                 | O(M)      |

Whole first; wrap and no jump exclude each other (wrap keeps the box filled, no jump keeps paths continuous).

## 2. The box: `PeriodicBox` in `@molstar/core`

`math/geometry/periodic-box.ts`: `{ basis: Mat3 (columns a, b, c), inverse: Mat3, origin?: Vec3 }` with
`minimumImage(out, d)` (fractional, round, back: exact for rectangular boxes and for bonded neighbours in any box
GROMACS accepts) and `fromCell(cell)`. `Frame` gains `box?: PeriodicBox`; `cell` stays (symmetry and the unit cell
representation read it) and is derived from the box. Why a new type: `Cell` rebuilds the basis in the standard frame
and has no origin ([analysis §3.2](analysis.md#3-gaps-found-on-the-way)). Ask first: it is a new public API (Q2).

## 3. Whole: the unwrap plan

Built once per topology (static data on the root model, `UnwrapPlan.get(model)`): `order` (atoms, parents first),
`parent` (`-1` for a root), `molecule` (component index) and per molecule `start`/`count` in `order`.

- **With topology bonds** (`IndexPairBonds`: PSF, PRMTOP, TOP, LAMMPS data): a BFS forest over them. Exact.
- **Without** (GRO, PDB, multi-model files): from the hierarchy, not distances, because the topology frame may itself
  be split ([analysis §1](analysis.md#1-the-symptom-and-why)). In a residue each atom's parent is the previous atom;
  a residue's first atom's parent is the previous residue's first atom when both are polymer residues of the same chain
  and consecutive in sequence; otherwise it is a root and starts a molecule. Holds while hops are under half a box,
  which the minimum-image convention guarantees for real systems; a wrong guess costs a molecule left split, never a
  wrong bond.

Per frame: `for i of order: if parent[i] >= 0: x[i] = x[p] + minimumImage(x[i] - x[p])` into fresh `Float32Array`s
(conformations are immutable once a model holds them, `model.md §3`).

## 4. Where it lives

- **Library** (`packages/model`, `coordinates/periodic.ts`): `UnwrapPlan`, `applyWhole`, `applyWrap`, and
  `PeriodicTrajectory.create(trajectory, params): Trajectory`, a lazy wrapper: `getFrameAtIndex(i)` gets the source
  frame and returns `{ ...frame, id, atomicConformation: fixed }`; frames without a box pass through. The representative
  is fixed too, and gets `atomicRanges` recomputed: polymer gaps are computed from the conformation at model creation
  (`basic/parser.ts`) and shared by frames, so a split topology frame would keep its false gaps. Bonds need nothing:
  topology bonds are right once molecules are whole, and distance bonds are computed per structure from the fixed
  representative. Covers every source, and headless and MVS get it without the plugin.
- **Plugin**: a transform `TrajectoryPeriodicity` (Trajectory → Trajectory, params `{ whole, wrap: off | box | selection,
  center: selection, noJump }`), inserted by the trajectory presets when the representative has a box. To check first:
  that the hierarchy manager and the trajectory animation accept a trajectory under a trajectory (else: a param group
  on `ModelFromTrajectory`, Q3).
- **Prism**: an operation `trajectory.periodicity` (aliases "fix periodic boundary", "make molecules whole", "center
  protein in box") next to the animation feature, and a toggle in the trajectory panel. **MVS**: a later spec extension
  of the trajectory node.

## 5. Tests

- Quick (`packages/core`): `PeriodicBox.minimumImage` for rectangular and triclinic boxes, against brute force over the
  27 neighbouring images.
- Quick (`packages/model`): a synthetic chain and waters placed across a face of a rectangular and a dodecahedral box:
  after whole every bonded distance is under 2 Å, the result is the same for the input shifted by any lattice vector,
  and whole is idempotent; wrap puts every molecule center within half a box of the center; a topology-bonds plan and a
  hierarchy plan agree. A two-frame triclinic GRO string through `PeriodicTrajectory`: no false polymer gaps.
- Quick (`packages/plugin/core`): the transform in a state tree, frame changes, toggling off restores the input.
- Heavy: `itWithin` bound for one frame of 1M atoms. By hand: AdK ([analysis §4](analysis.md#4-datasets)) in the
  viewer, before and after, at every frame.

## 6. Storybook: `Prism › Viewer › Simulations`

A new group for simulation output: MD trajectories with their box, and other simulations' frames. Its stories live in
`packages/prism/viewer/src/_stories/simulations/` (one file per subgroup: Storybook takes the title per file).

| Story                                 | What                                                                                                                                                                                                | Data                                                                                                     |
| ------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| **Periodic box › Split chain** (new)  | a peptide and waters crossing the faces of a triclinic box; args: off / whole / wrap / no jump, unit cell shown, frame; the default shows the fix, `off` shows the bug                                 | a synthetic multi-frame GRO in `data/examples/trajectory/periodic/` with its generator (no network; also a test fixture) |
| **Periodic box › AdK in water** (new) | the real case: protein torn in every frame, before and after                                                                                                                                       | MDAnalysis `adk_oplsaa.gro` + `.xtc` fetched at run time (raw.githubusercontent.com sends CORS `*`), as `Viewer/Data` fetches RCSB |
| **Trajectories › …** (done)  | one per coordinates format, the playback card on each, a `play` arg: PDB + DCD, PRMTOP + NetCDF, PDB + XTC, PDB + TRR, PDB + LAMMPS dump, multi-frame GRO (RST7 has no coordinates provider in v6) | PDB + DCD and PRMTOP + NetCDF the repository's; the rest mdtraj's test data (one alanine dipeptide run in every format, LGPL), fetched at run time |
| **Simularium › Frames, Playing**      | **moved** (done) from `Viewer/Particles` (`SimulariumFrames`, `SimulariumPlaying`): a simulation's frames, not cryo-ET picks                                                                                | `small.simularium`, as now                                                                               |

Staying where they are: `Viewer/Data`'s NMR ensemble and multi-model PDB (deposited ensembles; README decision 19 keeps
one story per data type there), `Viewer/Playback` (the playback card; it loads a DCD only to have frames), `Viewer/Dynamic
PDB` (an extension example), `Internals/Performance/Animation`. Nothing refers to the moved story ids but comments.

The Simularium move and the Trajectories stories are done (`_stories/simulations/`). The Periodic
box stories land with each phase: Split chain with phase 1 (its `wrap` and `no jump` args with phases 2 and 3), AdK
with phase 1.

## 7. Phases

0. **Readers** (small, independent; the `Viewer/Simulations` stories of §6 are done): GRO triclinic box (nine values); LAMMPS frames get a box (bounds, tilt factors) and
   use image flags or unwrapped coordinates when present; `Frame.box` from every reader (XTC, TRR keep their basis).
1. **Whole**: core box, plan, library wrapper, plugin transform, tests, the Periodic box stories (§6). The bug fix users see.
2. **Wrap** around the box or a selection (protein centered, solvent around it). A selection of several molecules
   wraps around the first one's center, iteratively; a complex split across the box needs GROMACS-style clustering:
   documented as a limit.
3. **No jump**: shifts per molecule per frame, computed lazily up to the highest frame asked and cached (frames are in
   memory for every reader; only molecule centers are needed per step).
4. **Exposure**: Prism operation and panel, the assistant skill text, MVS, docs. Later: periodic images (a P1 supercell
   from the box, through the existing symmetry structure).

Each phase: an entry in `.v6/changelogs/`, `model.md §3–4` and `formats.md` updated in the same commit.

## Open questions

1. **Default.** Whole on by default when frames have a box? Recommended: the unfixed view is never what one wants, it
   costs one O(N) pass per frame, and distances within a molecule do not change (distances between molecules can: a
   different image).
2. **`PeriodicBox`**, a new core type, or basis and origin added to `Cell`?
3. **Where the plugin applies it**: a Trajectory → Trajectory node (recommended: one place, composes) or a param of
   `ModelFromTrajectory`?
4. **No bonds**: the hierarchy heuristic (§3), or periodic-aware distance bonds as well (finds bonds across the face, so
   ligands with several residues are handled too; more work in `intra-compute.ts`)?
5. **Storybook group** (§6): `Viewer/Simulations`, with only Simularium moved into it, or also the ensembles from
   `Viewer/Data` and the trajectory in `Viewer/Playback`?
