# Analysis

## 1. The symptom, and why

A frame from a periodic run can hold a molecule split across the box: some atoms near one face, the rest near the
opposite face. Mol\* has no notion of periodicity outside crystal symmetry, so:

- **Bonds are fixed at the topology frame and reused.** `Unit.tryRemapBonds` (`packages/model/src/model/structure/structure/unit.ts`)
  keeps the old intra-unit bonds when the frames share one `IndexPairBonds` (topology from PSF, PRMTOP, TOP, LAMMPS)
  and, without index pairs, whenever `dynamicBonds` is off (the default) or the bonds are `canRemap` (dictionary-based
  and sequenced: proteins, nucleic acids, water). A bond whose atoms were imaged apart in a later frame is drawn the
  length of the box.
- **Polymer gaps are fixed at model creation.** `atomic-ranges.ts` decides "consecutive residues not connected" once,
  on the topology frame; later frames stretch the cartoon/trace across the box.
- **If the topology frame itself is split** (common: a GRO or PDB written by the engine), distance-based bonds miss the
  cross-face bonds, so the molecule is also broken into pieces even where nothing is stretched.

## 2. What Mol\* already has

- Every trajectory reader parses the box: `Frame.cell` (`coordinates/coordinates.ts`) from XTC and TRR (box vectors,
  `Cell.fromBasis`), DCD (three heuristics), NetCDF (lengths and angles). `trajectoryFromModelAndCoordinates`
  (`model.ts`) sets a P1 `ModelSymmetry` per frame from it, so the unit cell representation already follows the frames.
- `Coordinates.hasCell`; `Cell.toFractional`/`fromFractional`.
- Nothing reads the box to fix coordinates. No minimum-image helper in `@molstar/core` (`math/geometry/spacegroup` has
  cell and operators only).

## 3. Gaps found on the way

1. **GRO drops triclinic boxes.** `gro/parser.ts` keeps the first three of up to nine box values and `gro.ts` builds a
   rectangular cell (`PI/2` angles). AdK below is a rhombic dodecahedron (60°, 60°, 90°): its GRO cell is wrong, its
   XTC cell right.
2. **`Cell` loses the box's orientation and origin.** `Cell.fromBasis` keeps lengths and angles; `Cell.create` rebuilds
   the basis in the standard frame (a along x, b in xy). Right for GROMACS (lower-triangular boxes) and Amber; wrong for
   a DCD or LAMMPS box that is rotated or does not start at the origin. Making molecules whole needs only the lattice
   (orientation matters, origin does not); wrapping into "the box" needs the origin too.
3. **LAMMPS trajectory frames have no cell.** `lammps-trajectory.ts` uses `bounds` for scaled coordinates and the first
   frame's symmetry, but sets no `Frame.cell`, and ignores tilt factors (`xy xz yz`, triclinic). LAMMPS can also write
   image flags (`ix iy iz`) or unwrapped coordinates (`xu yu zu`), which are exact and better than any heuristic.
4. **DCD box origin is not in the file** (NAMD keeps it in `.xsc`); CHARMM/NAMD boxes are often centered on 0.

## 4. Datasets

Measured with Mol\*'s own XTC/GRO readers: per frame, residues whose atoms span more than 0.4 of the shortest box
vector, and the largest Cα–Cα distance between consecutive residues.

| Dataset                                                                                                                                                                           | Size                                 | Box                                | Result                                                                                                                                                              |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------ | ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **AdK in water**, MDAnalysis test data: [`adk_oplsaa.gro`](https://github.com/MDAnalysis/mdanalysis/raw/develop/testsuite/MDAnalysisTests/data/adk_oplsaa.gro) + `adk_oplsaa.xtc` | 2.1 MB + 1.7 MB; 47 681 atoms, 10 frames | rhombic dodecahedron, 80 × 80 × 57 Å | **split in every frame**: 5–11 protein residues split, consecutive Cα up to 79.5 Å apart (one box length); water whole. Also exercises the triclinic GRO bug (§3.1). |
| Alanine dipeptide, mdtraj test data: [`frame0.xtc`](https://github.com/mdtraj/mdtraj/raw/master/tests/data/frame0.xtc) + `native.pdb`                                              | 72 KB; 22 atoms, 501 frames          | triclinic                          | never split (stripped, already whole): a negative control and a "drifts out of the box" case for wrapping.                                                          |

AdK is the demo: small, public, triclinic, a protein visibly torn in each frame. The tests do not ship it (licence and
size; tests never reach the network): they use a synthetic box (plan §5). Licences of both test sets to be checked
before any fixture is derived from them.

## 5. How other tools do it

| Tool       | Operations                                                                                                                       |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------- |
| GROMACS    | `trjconv -pbc mol` (whole, centers in box), `whole`, `nojump`, `cluster`; `-center` on a group                                    |
| VMD        | `pbc unwrap` (no jumps over time), `pbc wrap -compound res/fragment -center com -centersel`, `pbc join` (whole)                    |
| MDAnalysis | on-the-fly transformations: `unwrap` (whole, by bonds), `wrap(compound=…)`, `center_in_box`, `nojump`                             |
| NGL        | trajectory params `removePbc` (whole, by atom order), `centerPbc`, `removePeriodicity` (no jumps between frames), `superpose`       |

Common ground: three separable operations. **Whole** (each molecule contiguous; per frame, needs bonds or an order of
atoms and the lattice). **Wrap** (each molecule's center into a box around a chosen center; per frame, needs whole
molecules). **No jump** (a molecule never jumps a box length between frames; needs the previous frame, so sequential).
