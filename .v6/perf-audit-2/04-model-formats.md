# Slice 04 v2: model formats and model

Scope: [`packages/model/src/formats`](../../packages/model/src/formats) (structure, volume, shape and topology
converters, format properties) and [`packages/model/src/model`](../../packages/model/src/model) without
`model/structure` (sequence, alignment, volume grid, loci, custom properties); the kinemage reader
([`packages/extensions/kinemage/src/reader`](../../packages/extensions/kinemage/src/reader)) is included because nothing else covers it.
HEAD `0944823d4`, 2026-10-22, Node 25.6, Apple M4, moderate shared load (±10–20 % on times). v1:
[../perf-audit/04-model-formats.md](../perf-audit/04-model-formats.md). Method and fixtures: [00-profile.md](00-profile.md). No
finding touches Canvas3D or `gl`. MVS loading is in [mvs-optimizations](../mvs-optimizations/research/README.md).

## Hot paths now

The main path is building a `Model` from parsed mmCIF/BCIF. It is cheap for small entries (7qpd: 4–6 ms) and large at
scale: flat-192 (1.12M atoms, BCIF) 0.48–0.58 s / 468–500 MB allocated; flat-864 (5M) 2.2 s / 2.1 GB.

| Path | Cost | Where it goes |
| --- | --- | --- |
| single-model BCIF, 1.12M atoms (flat-192) | 0.48–0.58 s | GC 13–22 %; `arrayView` (sorted `Table.view` copies, P5) 14–17 %; BCIF decode of every `atom_site` column (`getNormalizedAtomSite`) ~20 %; `_makeBuckets` 5 %; `addMissingResiduesFromEntityPolySeq` 3 % (P15); then hierarchy, index, derived data |
| sorted models, real entries | 9rxg 48 %, 21nr 79 % of rows moved | entities or chains not contiguous in the file, so the whole table goes through `Table.view` (P5), not only the fixture |
| BCIF model retained | ~320 B/atom (flat-192: 360 MB) | decoded source columns 188 B/atom + the model's sorted copies and conformation 133 B/atom (P16); text CIF 54 B/atom on top of its tokens |
| multi-model text mmCIF, 50 × 11.5k | 171–182 ms, 154–185 MB alloc | every hierarchy built, then compared (P4) |
| multi-frame GRO, 1000 × 327 atoms | 200 B/atom/frame, 65.4 MB | a full model per frame (P2) |
| IHM ensemble, 100 × 3.5k elements | ~720 KB/model retained | nothing shared between models (P3) |
| MOL2, 1000 poses × 26 atoms | 39–53 ms, 27–33 KB/pose retained | fixed per-model overhead per pose (P2) |

Already good: text tokenizer columns, hierarchy sharing of multi-model files (`sameAsPrevious`), particles, the
isosurface bounding sphere (memory order, cached), on-demand grid gradients, one Needleman-Wunsch in core (identical
chains aligned once), typed `IndexPairBonds`, the struct_conn partner lookup, IHM group names. The +30 % model
allocation the v2 profile shows on 7qpd does not reproduce in isolation (A/B against the v1 sources: 3.2–4.2 → 3.4–4.2
MB; flat-32 103 → 82 MB, 72 → 64 ms).

Changed since v1: #2017 fixed P1, P7–P12, P14 and the in-order half of P5; #2007 fixed P6 in `basic/properties.ts`;
#1989 fixed B2.

## v1 → now

27 entries (04/P1–P14, 04/B1–B7, and the 12/13 rows on this slice's files): **fixed 9, partly fixed 4, still open 13,
regressed 0, dropped 1**. Every v1 probe still open was rerun at HEAD and gives the v1 numbers.

| v1 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P1 | Needleman-Wunsch keeps three `(n+1)(m+1)` matrices | fixed | #2017 | `alignment.ts:9` wraps core `globalAlign` (two rows, 1-byte traceback); `pdb/atom-site.ts:169-182` aligns identical chains once |
| P2 | Multi-frame GRO builds a full model per frame | still-open | | `gro.ts:137`; rerun 200 B/atom/frame, 65.4 MB on multi1000.gro |
| P3 | IHM ensembles rebuild everything per model | still-open | | `parser.ts:280`; rerun +72.4 MB / 100 models, all released with the trajectory (checked with WeakRefs) |
| P4 | Multi-model mmCIF builds each hierarchy, then compares it | still-open | | `atomic.ts:237-241`; rerun 171–182 ms, 154–185 MB on multi50.cif |
| P5 | `sortAtomSite` copies all columns through `Table.view` | partly-fixed | #2017 | in-order ranges use windows (`sort.ts:47-60`); sorted models still `Table.view` everything (`:63`); `getNormalizedAtomSite` still decodes every column |
| P6 | Missing-residue scan builds a key per atom | partly-fixed | #2007 | `basic/properties.ts:83-99` keys per residue; `pdb/entity.ts:234` still per atom |
| P7 | `getIsosurfaceBoundingSphere` full scan with `getCoords` | fixed | #2017 | `grid.ts:97-164`, memory order, first/last hit per row, WeakMap cache |
| P8 | Grid gradients cached for the grid's lifetime | fixed | #2017 | `grid.ts:338` on demand from strides; no cache |
| P9 | `IndexPairBonds` edge props in holey arrays | fixed | #2017 | `index-pair.ts:44-50` typed arrays |
| P10 | struct_conn partner lookup tries every entity | fixed | #2017 | `struct_conn.ts:127-147` asym → entity ids |
| P11 | `guessElementSymbolString` regex per atom | fixed | #2017 | char codes for ASCII |
| P12 | `getModelGroupName` scans per model | fixed | #2017 | `basic/util.ts:20-37` two Maps once |
| P13 | Cube and DX volumes as Float64 | still-open | | `cube.ts:21,34`, `dx.ts:18` |
| P14 | `getResidueSourceIndex` sorts when in order | fixed | #2017 | `atomic.ts:121-132` increasing fast path |
| B1 | `getUniqueComponentNames` `memoize1` leak | still-open | | `properties.ts:189`; rerun +2.6–2.8 MB heap per load of multi2.cif |
| B2 | Dynamic format properties undefined for models 2..N | fixed | #1989 | `property.ts:76` presence in the store; assets per owner; `multi-model-props.test.ts` |
| B3 | ADP mapping drops the last atom / ids not 1..N | still-open | | `anisotropic.ts:53-63`; rerun `[0,1,-1]`, `[-1,-1,-1]` |
| B4 | `Sequence.ofSequenceRanges` off by one | still-open | | `sequence.ts:212`, `properties/sequence.ts:121`; rerun seqId(0) = 2, last = 21 |
| B5 | CONECT triple bonds written as `'trpl'` | still-open | | `conect.ts:131` |
| B6 | Granularity options under 'With Symmetry' | still-open | | `loci.ts:304`; rerun: residue..operator in 'With Symmetry' |
| B7 | `trajectoryFromModelAndCoordinates` frames drop `IndexPairBonds` | still-open | | `model.ts:121-128`; rerun: model 11 bonds, frames undefined |
| 13/B1 | Text mmCIF stops at ~4.8M atoms | dropped | | outside the slice (io); the app path decodes into a `ChunkedBigString` and `CIF.parse` takes `StringLike`; not re-run at 5M |
| 13/P3 | Model stage at 5M unshared | partly-fixed | #2007, #2017 | flat-864 2933 → 2221 ms, 2696 → 2120 MB; remainder is P5, P15, P16 |
| 12/B1 | Isosurface and segmentation bounds: `getCoords` per voxel | partly-fixed | #2017 | grid part fixed; `segmentation.ts:101-113` unchanged |
| 12/T1 | Trilinear interpolation through a `get()` closure | still-open | | `grid.ts:214`, 8 calls per sample |
| 12/L2 | Accessor call per value (CCP4 writer, cell bounds, periodic mapping) | still-open | | `ccp4.ts:99`, `volume.ts:359,372,562` |
| 12/B3 | Periodic trilinear interpolation reads unwrapped indices | still-open | | `grid.ts:250` |

## Performance findings (open, ranked)

| Id | vs v1 | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P5 | partly-fixed | Sorted models copy every `atom_site` column; normalization decodes every column | [sort.ts:63](../../packages/model/src/formats/structure/basic/sort.ts#L63) | time+memory | medium | flat-192 −45 ms of ~0.5 s, −90 MB alloc; flat-864 ~−0.2 s, −0.4 GB | S | confirmed | no |
| P16 | new | A BCIF model retains ~320 B/atom (string columns as pointer arrays, decoded cache, sorted copies) | [decoder.ts:219](../../packages/io/src/common/binary-cif/decoder.ts#L219) | memory | medium | −100–150 B/atom: −110–170 MB at 1.1M, −0.5–0.75 GB at 5M | L | confirmed | no |
| P4 | still-open | Multi-model mmCIF builds each hierarchy before comparing | [atomic.ts:241](../../packages/model/src/formats/structure/basic/atomic.ts#L241) | time | medium | ~180 → ~105 ms per 576k rows (text) | M | confirmed | no |
| P2 | still-open | Multi-frame GRO (and multi-molecule MOL2) build a full model per frame | [gro.ts:137](../../packages/model/src/formats/structure/gro.ts#L137) | memory | medium | 200 → ~15–25 B/atom/frame; 10k poses ~330 MB → a few MB | S | confirmed | no |
| P3 | still-open | IHM ensembles rebuild coarse hierarchy, sequence, ranges, static props per model | [parser.ts:280](../../packages/model/src/formats/structure/basic/parser.ts#L280) | time+memory | medium | ~720 → ~60–100 KB/model; ~2× faster load | M | confirmed | no |
| P15 | new | `missingResidues` computed eagerly on every load | [parser.ts:49](../../packages/model/src/formats/structure/basic/parser.ts#L49) | time | low | −18–31 ms at 1.1M, −150–190 ms at 5M | S | confirmed | no |
| P17 | new | Kinemage parser: a tag array per point, JS arrays for positions | [kinparser.ts:114](../../packages/extensions/kinemage/src/reader/kinparser.ts#L114) | memory | low | 316 → 208 B/point shared tags; ~60–80 with typed arrays | S | confirmed | no |
| P19 | new | `molstar_bond_site` bonds rebuilt per model over every model's rows, and per coordinate frame | [mmcif.ts:98](../../packages/model/src/formats/structure/mmcif.ts#L98) | time | low | O(M·(A + M·B)) → O(A + M·B) per file | S | plausible | no |
| P13 | still-open | Cube and DX volumes as Float64; cube orbitals copied again | [cube.ts:21](../../packages/model/src/formats/volume/cube.ts#L21) | memory | low | −50 % volume memory | S | confirmed | no |
| 12/T1 | still-open | Trilinear interpolation: 8 `get()` calls through a closure per sample | [grid.ts:214](../../packages/model/src/model/volume/grid.ts#L214) | time | low | ~1.4× per sample | S | plausible | no |
| P6 | partly-fixed | PDB SEQRES unobserved-residue scan: template key per atom | [entity.ts:234](../../packages/model/src/formats/structure/pdb/entity.ts#L234) | time | low | few ms per 100k-atom PDB | S | confirmed | no |
| P18 | new | Volume stats: four full passes (min, max, mean, rms) | [ccp4.ts:118](../../packages/model/src/formats/volume/ccp4.ts#L118) | time | low | 41–52 → 13–26 ms on 256³ | S | confirmed | no |
| 12/B1 | partly-fixed | Segmentation bounds: `getCoords` per non-zero voxel | [segmentation.ts:104](../../packages/model/src/formats/volume/segmentation.ts#L104) | time | low | ~3–10× for that pass | S | confirmed | no |
| 12/L2 | still-open | `Volume.Cell.getBoundingSphere` and CCP4 writer: accessor per cell | [volume.ts:359](../../packages/model/src/model/volume/volume.ts#L359) | time | low | ~2× per call | S | confirmed | no |

### P5: sorted `atom_site` copies

- Mechanism: a model whose rows moved goes through `Table.view(atom_site, atom_site._schema, indices)`; `arrayView`
  copies every array-backed column, including those the model never reads (`pdbx_label_index`, `pdbx_sifts_xref_db_*`
  ×4), and `Cartn_x/y/z` as Float64, which `getConformation` (`atomic.ts:150-152`) copies again to Float32 and drops.
  Before that, `getNormalizedAtomSite` (`basic/util.ts:55`) runs `Table.ofColumns` over the lazy `CategoryTable`, which
  creates every field, and a BCIF field decodes when it is created (`binary/field.ts:15`).
- Evidence: `view.ts` on flat-192 (rerun): all 26 present columns 88–92 ms / 205 MB; only the 21 used 64–74 ms /
  170 MB. Unused-column decode 22 ms of 147 ms (`cols.ts`). Not fixture-only: 9rxg 48 %, 21nr 79 % of rows move
  (`sorted.ts`). The unused copies are garbage after the model is built (alloc and time, not retained).
- Fix: (1) the sorted branch builds one memoized lazy getter per column (copy on first read); (2) `getConformation`
  gathers Cartn straight to Float32 through the source index; (3) `getNormalizedAtomSite` replaces only the five
  substituted columns over the original getters.
- Risk: hierarchy fast paths need plain-array columns (#2017's note) and `isHierarchyDataEqual` compares them; pin with
  `sort.test.ts`.

### P16: BCIF retained memory

- Mechanism: `stringArray` decodes to a JS `string[]` per column (8 B/atom in Node, 4 in Chrome); the category caches
  every decoded field for the life of the source data (`binary/parser.ts:38`, `cif/schema.ts:208`). A sorted model then
  keeps `arrayView` copies of atom/comp ids, alt id and type symbol (6 × 8 B), Int32 charge and id, Float64 occupancy and
  B. Coordinates live twice: Float64 in the source, Float32 in the conformation.
- Evidence: `retain2.ts` rerun, flat-192: source columns 188 B/atom, the model on top 133 B/atom (321 B/atom = 360 MB
  at 1.12M; text CIF 54 B/atom). The saving of the fix is estimated, not measured.
- Fix: keep BCIF string columns dictionary-encoded (`value(i) = strings[indices[i] + 1]`, `areValuesEqual` on
  indices) so sorting gathers 4-byte indices; occupancy and B as Float32 in the conformation. Needs a design note
  (io + core `Column`).
- Risk: readers of `column.__array` as strings; `sort.ts:56`'s typed-array test; `Table.areEqual` fast paths.

### P4: multi-model hierarchy comparison

- Mechanism: `getAtomicHierarchy` runs `findHierarchyOffsets`, `createHierarchyData` and
  `createChainOperatorMappingAndSubstituteNames` (a `toDatabase` of `molstar_atom_site_operator_mapping` per model)
  before `isHierarchyDataEqual`; text CIF compares go through `str`/`int` and allocate. NMR/MD files throw it away.
- Evidence: `multi.ts` rerun, multi50.cif: 171–182 ms, 154–185 MB alloc (v1 176–189 ms; the compare was 65 ms).
- Fix: compare offsets with the previous segments, then the hierarchy source columns row by row on `data.atom_site`
  (`areValuesEqual(prevRow, curRow)`, raw-token for text); build only on a mismatch. Read the operator mapping once per
  file.
- Risk: the compared set must equal what `createHierarchyData` builds; test multi-model PDB, BCIF, differing counts.

### P2: GRO frames, MOL2 poses

- Mechanism: `trajectoryFromGRO` calls `getBasic` + `createModels` per frame with no previous model (TODO at `:129`);
  `mol2.ts` does the same per molecule.
- Evidence: `grobench2.ts` rerun: 200 B/atom/frame (coordinates alone 12); `mol2many.ts`: 1000 poses 27–33 KB each.
- Fix: one `BasicData` for all frames (`pdbx_PDB_model_num` = frame) and one `createModels`, then `ModelSymmetry` per
  frame; MOL2 the same for consecutive molecules with identical atom tables.
- Risk: differing topologies keep their own hierarchy (the compare handles it); per-pose charges and bonds stay per
  model. Sharing a hierarchy hits B9.

### P3: IHM ensembles

- Mechanism: `createIntegrativeModel` per `ihm_model_list` row: `getCoarse`, `getSequence` over all
  `entity_poly_seq`, `getAtomicRanges`, `getStructAsymMap`, fresh `customProperties` and `_staticPropertyData`; one
  `Mat3` per gaussian (`coarse.ts:89`).
- Evidence: `ihmmem2.ts` rerun: +72.4 MB for 100 models. v1's script read "leaked" after dropping the trajectory; with
  an event-loop turn before GC (`verify/ihmref.ts`) all of it is freed, so it is per-model data, not a cache.
- Fix: compare sphere/gaussian site keys and the atomic part with the previous model; reuse hierarchy, sequence,
  ranges, entities, properties, static data; `getSequence` once per file; covariances as one `Float32Array(9n)`.
- Risk: the flat covariance changes `CoarseGaussianConformation` (graphics, export).

### Low

- P15: `getCommonProperties` calls `getMissingResidues` eagerly. Correction to the audit: it is read, by the classic
  plugin UI sequence panel ([polymer.ts:98](../../packages/plugin/ui/src/sequence/polymer.ts#L98), greys missing
  residues), so it cannot be dropped; Prism and Seq\* do not read it. `missing.ts` rerun: 18–19 ms on flat-192,
  147–188 ms on flat-864. Fix: make the `MissingResidues` object itself lazy (compute on first `has`/`get`/`size`); not
  a getter on `properties`, which `createStandardModel` spreads (a getter would be evaluated by the spread).
- P17: `pointMasters: string[] = []` per point grows to 17 slots on the first push; trimmed line, split array,
  position array and label per point; positions spread-pushed into `number[]`. `kin.ts` rerun: 500k dots + 500k
  vectors 664–666 ms, 316 B/point; shared tag arrays (`kinpatch`) 208 B/point, same time. Dead `vector3Positions` loop
  in `convertKinTriangleArrays` (`:262-265`, also wrong: `convertedPositions[i * 3] + 1`). Fix: cache tag arrays by
  tag string, typed/chunked accumulation, delete the dead loop.
- P19: `indexPairBondsFromMolstarBondSite` runs per model (`asDynamic`): `toDatabase`, a Map over the model's atoms,
  every `molstar_bond_site` row including the other models'. Per coordinate frame of an mmCIF topology it recomputes
  (and each visited frame keeps its own copy). By code; multi-model files with `molstar_bond_site` are rare. Fix: split
  rows per model once (WeakMap on the format), copy the topology model's value to frames; combine with B8.
- P13: `Tensor.Space(..., Float64Array)` in cube/dx; orbital extraction allocates another Float64. Fix: parse to
  Float32 in `io/reader/cube`, `io/reader/dx`.
- 12/T1: `getValue` closure, 8 calls with three-index math per sample (volume-value themes per vertex). Fix: base
  offset from `Tensor.strides` plus 7 deltas, as `makeGetGradient`; wrap when periodic (12/B3). Not re-measured.
- P6: `${modelN}|${labelAsym}|${labelSeq}` per polymer atom row in the PDB SEQRES path. Fix: emit on change only.
- P18: four passes in ccp4/dx/cube/dsn6. `stats.ts` on 256³: 41–52 ms vs one fused pass 13–26 ms. Fix: one
  `arrayStats` helper in core (also B10).
- 12/B1: segcif bounds call `getCoords` per non-zero voxel. Fix: memory-order loops (as `grid.ts:118-139`), with B13's
  guard.
- 12/L2: `getBoundingSphere` of a cell loci transforms each cell twice per instance in two passes; the CCP4 writer
  calls `dataOffset` per voxel. Fix: one transform per instance, cache positions between passes.

## Bugs

| Id | vs v1 | Severity | Title | File | Verdict |
| --- | --- | --- | --- | --- | --- |
| B8 | new | medium | `molstar_bond_site` bonds built in element order, read as source rows; models 2..N with their own topology get none | [mmcif.ts:106](../../packages/model/src/formats/structure/mmcif.ts#L106) | confirmed |
| B10 | new | medium | Volume sigma is RMS, not standard deviation (CCP4 without ARMS, DX, Cube, DSN6) | [ccp4.ts:122](../../packages/model/src/formats/volume/ccp4.ts#L122) | confirmed |
| B1 | still-open | medium | Leak: `getUniqueComponentNames` `memoize1` keeps every `BasicData` | [properties.ts:189](../../packages/model/src/formats/structure/basic/properties.ts#L189) | confirmed |
| B7 | still-open; **fixed** #2032 | medium | Coordinate-trajectory frames lose bonds set with `IndexPairBonds.Provider.set` | [model.ts:127](../../packages/model/src/model/structure/model/model.ts#L127) | confirmed |
| B3 | still-open; **fixed** #2033 | medium | ADP mapping drops the last atom and ids not 1..N; shared-hierarchy models get model 1's ADPs | [anisotropic.ts:53](../../packages/model/src/formats/structure/property/anisotropic.ts#L53) | confirmed |
| 12/B3 | still-open; **fixed** #2032 | medium | Periodic trilinear interpolation reads unwrapped indices | [grid.ts:250](../../packages/model/src/model/volume/grid.ts#L250) | confirmed |
| B9 | new | low | Models sharing the previous hierarchy also get its `atomSourceIndex` | [atomic.ts:244](../../packages/model/src/formats/structure/basic/atomic.ts#L244) | confirmed |
| B11 | new; **fixed** #2033 | low | One-sided PDB CONECT (higher serial only) dropped | [conect.ts:70](../../packages/model/src/formats/structure/pdb/conect.ts#L70) | confirmed |
| B4 | still-open | low | `Sequence.ofSequenceRanges` off by one; last coarse element dropped | [sequence.ts:212](../../packages/model/src/model/sequence/sequence.ts#L212) | confirmed |
| B5 | still-open; **fixed** #2033 | low | CONECT triple bonds written as `'trpl'` | [conect.ts:131](../../packages/model/src/formats/structure/pdb/conect.ts#L131) | confirmed |
| B6 | still-open; **fixed** #2033 | low | Granularity options grouped under 'With Symmetry' | [loci.ts:304](../../packages/model/src/model/loci.ts#L304) | confirmed |
| B12 | new | low | MOL2 bonds assume atom ids 1..N in file order | [mol2.ts:163](../../packages/model/src/formats/structure/mol2.ts#L163) | plausible |
| B13 | new | low | Segmentation load throws on a voxel label without a table set | [segmentation.ts:105](../../packages/model/src/formats/volume/segmentation.ts#L105) | plausible |
| B14 | new | low | Kinemage: a re-opened `@group` keeps the previous group's masters | [kinparser.ts:767](../../packages/extensions/kinemage/src/reader/kinparser.ts#L767) | plausible |

### B8: `molstar_bond_site` index space

- Mechanism: `idToIndex.set(atomId.value(i), i)` maps ids to element indices of the sorted conformation, but every
  `IndexPairBonds` reader treats indices as source rows (`intra-compute.ts:102` `offset[sourceIndex.value(aI)]`,
  `:453/:465` and `inter-compute.ts:160` `invertedIndex[b[i]]`), as `_spec/formats.md` §2 says. They differ whenever
  `sortAtomSite` reordered rows. For a later model with its own hierarchy, `atomSourceIndex` holds absolute file rows
  while `invertIndex` (`model.ts:215-222`) sizes by atom count: writes past the end, inverted index all zeros
  (`packages/model/TODO.md` already lists that half).
- Failure (rerun): `srcidx2.ts`, chains A and B interleaved (source index `[0,2,1,3]`), bonds 1–3, 2–4: 0 intra and 0
  inter bonds instead of 2. `srcidx.ts`, model 2 with 4 atoms (source index `[3,4,5,6]`): IndexPairBonds 3 edges, unit
  bonds 0. Reachable from Mol*'s own export: an assembly with ligands interleaves entities per operator, so reloading
  sorts it.
- Fix: map ids to model-local source rows (`atomSourceIndex.value(i) − modelStart`); settle one space for
  `atomSourceIndex` (or store the offset with the inverted index). Add sorted and changing-topology cases to
  `multi-model-props.test.ts`.

### B10: sigma as RMS

- Mechanism: `IsoValue` uses `sigma` as a standard deviation (`absolute = relative · sigma + mean`,
  `volume.ts:80-130`). The fallbacks use core `arrayRms` = `sqrt(sd² + mean²)`: `ccp4.ts:122` (ARMS NaN or ≤ 0, e.g.
  IMOD's −1), `dx.ts:34`, `cube.ts:61`, `dsn6.ts:53`. MTZ (`mtz.ts:202`) and structure factors (`:89`) use the sd; the
  CCP4 header's ARMS is itself the deviation from the mean.
- Failure (`stats.ts` rerun, 256³, mean 103, sd 7.35, max 117): sigma 103.26; 1σ = 206 > max, so the default isosurface
  is empty and the relative slider spans ~±0.1. Maps with mean ≈ 0 (most X-ray, APBS DX) are barely affected.
- Fix: one-pass `arrayStats` (min, max, mean, sd) in core for the four readers (also P18). Saved absolute isovalues of
  affected maps change (they were wrong).

### B1, B7, B3, 12/B3

- B1: module-level `Map` in `memoize1` keyed by `BasicData`, reached when `chem_comp` is empty; the key holds
  `atom_site` and through it the source text. `leak2.ts` rerun: +2.6–2.8 MB heap per load, never released. Fix: compute
  the names once in `getCommonProperties` and pass them, or a `WeakMap`.
- B7: frames get a fresh `_dynamicPropertyData` and only `ModelSymmetry` again; SDF/MOL/MOL2/CIF-core/LAMMPS bonds set
  with `Provider.set` have no registry entry. `trajbonds.ts` rerun: model 11 bonds, frames 0–2 undefined (distance
  bonding, orders lost). Fix: copy the model's `IndexPairBonds` onto each frame as `trajectoryFromTopologyAndCoordinates`
  does (`model.ts:166-174`). Adjacent slice (`model/structure`), listed in `packages/model/TODO.md`.
- B3: `new Int32Array(atomId.rowCount)` indexed by id value. `aniso.ts` rerun: ids 1..3 → `[0,1,-1]`, ids 101..103 →
  `[-1,-1,-1]`; the last atom never has an ellipsoid, PDB ANISOU files with TER gaps lose trailing atoms. New: the
  provider is static, and `sameAsPrevious` models share `_staticPropertyData`, so model k shows model 1's mapping. Fix:
  `Map<id, element>`, `ei !== undefined`; key the mapping per conformation.
- 12/B3: `isPeriodic` skips the bounds check but `i, j, k` are passed unwrapped and `ii, jj, kk` clamped with
  `Math.min`. Volume-value color/size themes on MTZ, structure-factor and whole-cell CCP4 maps: atoms outside the unit
  cell (common) read the wrong voxel or `undefined` → NaN → default color. Fix: wrap each index modulo its dimension.

### Low bugs

- B9: when `isHierarchyDataEqual` holds, `getAtomicHierarchy` returns `previous.atomicHierarchy`, whose
  `atomSourceIndex` is model 1's rows. `srcidx.ts`: model 2's source index `[0,1,2]`, its rows are 3..5. mmCIF export of
  model 2 writes `pdbx_label_index` and `pdbx_sifts_*` from model 1's rows (`atom_site.ts:179`); the label's `[idx N]`,
  MolScript `atom.sourceIndex` and SIFTS mapping (`sifts-mapping.ts:90`) read model 1's rows. Local consumers (inverted
  index, frame coordinates) work by accident. Fix together with B8: per-model offset, not a reinterpretation.
- B11: `if (idxA > idxB) continue` drops a bond listed only on the higher serial's CONECT. `conect.ts` probe rerun:
  `CONECT 2 1` alone → 0 struct_conn entries. Fix: dedupe pairs with `sortedCantorPairing`, count repeats per pair, take
  the max multiplicity.
- B4: `row + minSeqId + 1` (should be `row + minSeqId`) and `offsets[cI + 1] − 1` as an exclusive end.
  `seqranges.ts` rerun: seqId(0) = 2, last = 21. Coarse models without `entity_poly_seq` only.
- B5: `'trpl'` instead of `'trip'`; `CONECT 1 2 2 2` reads back as order 1.
- B6: `k.indexOf('Instances') ? …` is truthy for −1; use `includes`.
- B12: `indexA/B = origin/target_atom_id − 1`, used as rows. A MOL2 numbered 1,2,4,5 bonds the wrong pairs or indexes
  past the atom count. Plausible (MOL2 writers usually number 1..N); no fixture. Fix: id → row Map.
- B13: `setBounds[v]` exists only for set ids in `segmentation_data_table` with a non-zero segment; any other voxel
  label throws (`b[0]` of undefined) and the volume fails to load. Plausible; no fixture. Fix: `if (!b) continue`.
- B14: `currentGroupMasters = groupMasters` sits inside `if (!kinemage.groupDict[groupName])`, so a second `@group {A}`
  keeps the previous group's masters for its lists. Plausible; no fixture. Fix: keep a group's masters in `groupDict` and
  restore them on re-open.

## Checked and dropped

Nothing refuted. Corrections:

- P15: "read by nothing in the repo" is wrong; the classic UI sequence panel reads it. Make it lazy, do not remove it.
- P3: v1's `ihmmem2.ts` printed "+72.4 MB leaked after dropping the trajectory"; that was measured before the event loop
  turned. With WeakRefs and one turn the trajectory and all 100 models are collected (+0.8 MB left), so the 720 KB/model
  is per-model data, as v1 concluded.
- P5: the unused-column copies are freed after the build (alloc and time); the retained part belongs to P16.
- 13/B1: outside this slice (io); not re-verified at 5M.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/model-formats/`,
bundled with `node build.mjs <x.ts> <x.mjs>` (esbuild, conditions `molstar-src`, `import`, `default`), run with
`node --expose-gc x.mjs <file>`. Fixtures: `tests/build/pipeline/fixtures/1tup-flat-{32,192,864}.bcif`; v1 inputs
(`multi50.cif`, `multi2.cif`, `multi1000.gro`, `ihm100.cif`) in the v1 scratchpad
`/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/model-formats/`.

| Script | Finding | Input |
| --- | --- | --- |
| `view.ts`, `cols.ts`, `sorted.ts`, `modelprof.ts`, `sumprof.mjs`, `alloc.ts` (+ `buildv1.mjs`) | P5, hot paths | flat-192, flat-32, 9rxg, 21nr, 7qpd |
| `retain.ts`, `retain2.ts` | P16 | flat-192, flat-32 |
| `missing.ts` | P15 | flat-192, flat-864 |
| `mol2many.ts` | P2 | synthetic 1000 poses |
| `kin.ts`, `kin2.ts`, `kinpatch/` | P17 | synthetic 500k dots + 500k vectors |
| `stats.ts` | B10, P18 | synthetic 256³ |
| `srcidx.ts`, `srcidx2.ts` | B8, B9 | synthetic mmCIF |
| `conect.ts` | B11 | synthetic PDB |
| `seqres.ts` | P1 | `seqres.pdb` |
| `v1re/*` (`multi`, `multib`, `grobench2`, `ihmmem2`, `leak2`, `aniso`, `seqranges`, `coarseseq`, `gran`, `trajbonds`) | P2–P4, B1, B3, B4, B6, B7 | v1 inputs; `data/examples/long_animation.sdf` |
| `verify/ihmref.ts` | P3 | `ihm100.cif` |
