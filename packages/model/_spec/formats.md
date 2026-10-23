# Formats

`src/formats/`: readers into models, coordinates, topologies, volumes and particle lists; writers.

## 1. Structures

- `ModelFormat<T> = { kind, name, data }`; each format has a namespace with `is()` (by `kind`). PDB models are
  `MmcifFormat` with the PDB kept in `data.source`.
- Every structure reader builds `BasicData` (`basic/schema.ts`, a subset of mmCIF tables; `createBasic` fills missing
  tables and, normalizing, substitutes `label_*`/`auth_*` columns for each other) and calls
  `createModels(basic, format, ctx)` (`basic/parser.ts`), which returns an `ArrayTrajectory`:
  - standard: `atom_site` split by `pdbx_PDB_model_num`, each run sorted (`sortAtomSite`) and made a model; a model
    whose tables equal the previous one's shares its topology (model.md §3);
  - integrative (`ihm_model_list` present): atomic plus coarse spheres and gaussians.
- Readers: mmCIF/BinaryCIF (`mmcif.ts`; `MmcifFormat` sits apart in `mmcif-format.ts` so code can recognize mmCIF models
  without loading the parser), CCD (ideal and model coordinates), PDB/PQR/PDBQT (via `pdbToMmCif`), GRO, MOL, SDF, MOL2,
  XYZ, CIF core, cube, LAMMPS data and trajectory.
- Coordinates only: DCD, XTC, TRR, NetCDF, LAMMPS trajectory. Topologies only: PSF, PRMTOP, TOP (bonds as 0-based index
  pairs).

## 2. Format properties

`FormatPropertyProvider.create(descriptor, { asDynamic?, cachePerFormat? })` (`structure/common/property.ts`): a
registry from `ModelFormat.kind` to an `obtain(model)`.

- `get` stores the result, `undefined` included, in the static or dynamic bag under the descriptor's name and checks
  that bag, not `customProperties`: same-topology models share `customProperties` but not dynamic data, so the
  descriptor's presence would hide the value of every model after the first one asked.
- `cachePerFormat` memoizes per `model.sourceData` (a `WeakMap`): one value for every model of a file.
- Registration is an import side effect of the reader module (`mmcif.ts`, `cif-core.ts`; `pdb.ts` imports `mmcif.ts` for
  it).

| Provider                  | Bag                 | From                                                                         |
| ------------------------- | ------------------- | ---------------------------------------------------------------------------- |
| `ModelSymmetry`           | dynamic, per format | mmCIF (assemblies, spacegroup, NCS); `fromCell` by other readers, frames     |
| `Model.PeriodicBox`       | dynamic             | frames with a cell or `Frame.box`; GRO; LAMMPS; NAMD `.xsc` (`namd-xsc.ts`)  |
| `ModelSecondaryStructure` | static              | mmCIF `struct_conf`, `struct_sheet_range`                                    |
| `StructConn`              | static              | mmCIF `struct_conn`                                                          |
| `ComponentBond`           | static, per format  | mmCIF `chem_comp_bond`; set by CCD                                           |
| `IndexPairBonds`          | dynamic             | mmCIF `molstar_bond_site`; set by MOL/SDF/MOL2, CIF core, LAMMPS, topologies |
| `AtomSiteAnisotrop`       | dynamic             | mmCIF (by each model's atom ids), CIF core (by atom label)                   |
| `AtomPartialCharge`       | static              | set by PDBQT, MOL2, LAMMPS                                                   |

`IndexPairBonds` index atoms by source index (`atomSourceIndex`); `maxDistance` (default -1: an element-based threshold;
`Infinity`: trust the table; a per-bond `distance` takes precedence), `cacheable`, and operators per bond (assemblies
written with operator ids). The frames of `trajectoryFromModelAndCoordinates` share the model's value: no format
provider gives a value a reader set.

## 3. Volumes and particles

- Volumes: CCP4/MRC, DSN6, DX, cube, density server, segmentation, structure factors and MTZ (through
  `shared/reflections.ts`). Each gives a `Volume` with one identity instance. CCP4 header stats are trusted unless they
  are all zero or NaN.
- Particles: Artiatomi, Dynamo, RELION STAR, CryoET Data Portal NDJSON, Simularium (a lazy `ParticleTrajectory`), mmCIF
  assemblies (`createParticleListFromMmcifAssembly`).

## 4. Writers and export

- `writer/`: MOL, SDF and MOL2 encoders over `atom_site` (component atoms and bonds supplied by the caller), and
  `CCP4Writer.writeMrc`.
- `model/structure/export/mmcif.ts`: `to_mmCIF`, `encode_mmCIF_categories` for structures of one model. Categories are
  copied from the source or generated (`_entity`, `_atom_site`, secondary structure, `_struct_asym`, `_entity_poly*`,
  `_chem_comp*`, `_pdbx_nonpoly_scheme`), plus custom properties with `cifExport`. Mol\* extensions:
  `molstar_atom_site_operator_mapping` (read back by the parser) and, if asked, `molstar_bond_site`.
