# Slice 04: model formats and model

Scope: [`packages/model/src/formats`](../../packages/model/src/formats) (structure, volume, shape and topology
converters, format properties) and [`packages/model/src/model`](../../packages/model/src/model) without
`model/structure` (sequence, alignment, volume grid, loci, custom properties). Fixtures and method:
[00-profile.md](00-profile.md). Measured 2026-10-19, Node 25.6, Apple M4. No finding touches Canvas3D or `gl`. MVS
loading is in [mvs-optimizations](../mvs-optimizations/research/README.md).

## Hot paths

On the profiled path (one model of mmCIF/BCIF) the slice is cheap: the `model` stage is 7 ms for 7qpd (11.5k atoms) and
90 ms for big-bcif (186k atoms), under 1 % of the pipeline.

| Path | Cost | Where it goes |
| --- | --- | --- |
| single-model BCIF, 186k atoms | 90 ms | decode of every `atom_site` column forced by `getNormalizedAtomSite` (20 %), `Table.view` copies in `sortAtomSite` (`arrayView` 15 %), per-atom keys in `addMissingResiduesFromEntityPolySeq` (10 %), atomic index / derived data (other slice) |
| multi-model text mmCIF, 50 × 11.5k | 176–189 ms | `getAtomicHierarchy` 99 ms, of which `isHierarchyDataEqual` 65 ms (`str10` 35 ms); `getConformation` 28 ms |
| PDB with SEQRES | +9 ms, +40 MB alloc for one 555-residue chain | Needleman-Wunsch on `number[][]` (P1); per chain, scales with residues² |
| multi-frame GRO | 100–200 B/atom/frame retained | a full model per frame (P2) |
| IHM ensemble, 100 × 3.5k elements | 225–247 ms, ~720 KB/model retained | nothing shared between models (P3) |
| CPU isosurface / wireframe / dots | +35–120 ms per rebuild on 256³ | full-grid scan for the bounding sphere (P7); the GPU path skips it |

Already good: text tokenizer columns (zero-copy views, raw-token `areValuesEqual`), hierarchy sharing between models of
multi-model mmCIF (`sameAsPrevious`), particles (typed arrays, chunked updates), LAMMPS / XYZ trajectories (shared
topology).

## Performance findings (confirmed and plausible, ranked)

| Id | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen | Addressed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | Needleman-Wunsch keeps three `(n+1)(m+1)` `number[][]` matrices | [alignment.ts:32](../../packages/model/src/model/sequence/alignment/alignment.ts#L32) | time+memory | high | 2.3–5× faster, 5× less memory (flat Int32); ~50× less with 2-row + 2-bit traceback | M | confirmed | no | #2017 |
| P2 | Multi-frame GRO builds a full model per frame | [gro.ts:134](../../packages/model/src/formats/structure/gro.ts#L134) | memory | high | 100–200 → ~15–25 B/atom/frame (8–15×); faster load | S | confirmed | no |  |
| P4 | Multi-model mmCIF builds each hierarchy, then compares it through allocating string compares | [atomic.ts:217](../../packages/model/src/formats/structure/basic/atomic.ts#L217) | time | medium | multi-model text mmCIF model stage ~180 → ~105 ms per 576k atom rows | M | confirmed | no |  |
| P3 | IHM ensembles recompute and keep coarse hierarchy, sequence, ranges and static props per model | [parser.ts:227](../../packages/model/src/formats/structure/basic/parser.ts#L227) | time+memory | medium | ~720 → ~60–100 KB/model for identical topologies; load ~2× faster | M | confirmed | no |  |
| P8 | Grid gradients, 12 B/voxel, cached on the Grid for its lifetime | [grid.ts:310](../../packages/model/src/model/volume/grid.ts#L310) | memory | medium | −12 B/voxel retained once streamlines run (1.6 GB at 512³) | S | confirmed | no | #2017 |
| P7 | `getIsosurfaceBoundingSphere` scans the full grid with `getCoords` per hit | [grid.ts:85](../../packages/model/src/model/volume/grid.ts#L85) | time | low | 3–10× faster per call (35–120 → 11–17 ms on 256³); CPU mesh / wireframe / dots only | S | confirmed | no | #2017 |
| P5 | `sortAtomSite` copies all columns through `Table.view`, even for a contiguous identity range | [sort.ts:47](../../packages/model/src/formats/structure/basic/sort.ts#L47) | time | low | −13 ms of 90 (big-bcif), −24 ms of 137 (50-model BCIF); ~0.15 MB/model retained, not MBs | S | confirmed | no | #2017 |
| P6 | `addMissingResiduesFromEntityPolySeq` builds a template-string key per atom | [properties.ts:83](../../packages/model/src/formats/structure/basic/properties.ts#L83) | time | low | −8–12 ms on 186k atoms | S | confirmed | no | #2007 |
| P11 | `guessElementSymbolString`: regex and `toUpperCase` per atom | [util.ts:65](../../packages/model/src/formats/structure/util.ts#L65) | time | low | 100 → 26 ms per 2M atoms | S | confirmed | no | #2017 |
| P9 | `IndexPairBonds` graph stores operators, distance and flag in holey JS arrays | [index-pair.ts:44](../../packages/model/src/formats/structure/property/bonds/index-pair.ts#L44) | memory | low | −36 MB per 1M bonds in Node; −12–16 MB in Chrome (pointer compression) | S | confirmed | no | #2017 |
| P14 | `getResidueSourceIndex` comparator-sorts residues even when the order is identity | [atomic.ts:121](../../packages/model/src/formats/structure/basic/atomic.ts#L121) | time | low | ~1 ms of 5.4 (1tup-2x2x2); O(n) when monotonic | S | confirmed | no | #2017 |
| P13 | Cube and DX volumes stored as Float64; each cube orbital copied again as Float64 | [cube.ts:21](../../packages/model/src/formats/volume/cube.ts#L21) | memory | low | −50 % volume memory, needs the io readers to parse to Float32 | S | confirmed | no |  |
| P10 | struct_conn partner lookup tries every entity | [struct_conn.ts:139](../../packages/model/src/formats/structure/property/bonds/struct_conn.ts#L139) | time | low | rows × 2 × entities Map misses → rows × 2; ~10–20 ms on 10k-row / 100-entity entries, once per hierarchy | S | plausible | no | #2017 |
| P12 | `getModelGroupName` scans `ihm_model_group_link` per model | [util.ts:11](../../packages/model/src/formats/structure/basic/util.ts#L11) | time | low | O(models²) → O(models); matters above ~1k models | S | confirmed | no | #2017 |

### P1: sequence alignment matrices

- Mechanism: `S`, `V`, `H` are arrays of JS arrays grown cell by cell (`initMatrices`, `:49`); `-Infinity` makes rows
  double arrays; `scoreFn` closure with a string compare per cell; `trace` reads the three matrices. Callers:
  `alignCompIdsToSeqres` ([atom-site.ts:110](../../packages/model/src/formats/structure/pdb/atom-site.ts#L110)), once
  per chain on every PDB file with SEQRES (`computeSeqresAlignments`, `:169`); `AlignSequences.compute`
  (sequence superposition).
- Evidence (`align.ts`, `align2.ts`, rerun): n = 1000 / 2000 / 3000 → 24–28 / 91–114 / 204–207 ms and 50 / 224 /
  501 MB of heap. Flat `Int32Array` with the same recurrences and identical scores: 10–15 / 17–22 / 37–41 ms, 11 / 45 /
  101 MB. End to end (`pdbbench.ts`, `seqres.pdb`, 9010 atoms, one 555-residue chain): `trajectoryFromPDB` 7–8 ms /
  5–8 MB alloc without SEQRES, 16–23 ms / 46–49 MB with it. Cost is per chain and quadratic in chain length.
- Fix: flat `Int32Array` S/V/H (stride `m + 1`, `-0x3fffffff` sentinel), or two rolling rows plus a `Uint8Array`
  traceback with 2 bits per matrix; intern residue names to small ints once (or index a substitution matrix precomputed
  into an `Int8Array`); in `computeSeqresAlignments` memoize by (SEQRES array, observed key) so identical chains align
  once; fast path when observed is an in-order exact subsequence of SEQRES (greedy O(n + m)).
- Risk: `trace` tie-breaking (S, then V, then H) must be kept or ambiguous alignments change; compare `aliA` / `aliB`
  with the old implementation on random inputs and run the alignment and PDB `label_seq_id` tests.
- Addressed (#2017): one DP for the model and Seq\*, `@molstar/core/data/util/needleman-wunsch` (`globalAlign`): two
  Float64 rows, a one-byte traceback per cell (2 bits per state), a branch-free tie choice (mispredicted branches were
  ~40 % of the time). Both callers bit-identical (160k seeded cases against the old code; the model's leading gaps
  (`flat`) and extend-first ties are options). 3000 res: 210 ms / 1017 MB peak → 48 ms / 32 MB (4.4×, 32×; the estimate
  was 41 ms / 101 MB for flat Int32 matrices); 10000: 3.1 s / 3 GB → 0.5 s / 97 MB; Seq\* 3000: 82 → 45 ms.
  `seqres.pdb`: `trajectoryFromPDB` 17–29 ms / 45–55 MB → 10–21 ms / 7–19 MB. Identical chains align once; the greedy
  subsequence fast path was not done (it could change ambiguous alignments). Interning names to codes is done.

### P2: GRO trajectories

- Mechanism: `trajectoryFromGRO` calls `getBasic` + `createModels` per `gro.structures[i]` and never passes a previous
  model, so `sameAsPrevious` never applies: per frame, atoms / residues / chains tables, segmentations, atomic index,
  derived data, ranges, sequence, chemical component map, entities and a `guessElementSymbolString` pass. The TODO at
  `:126` says so.
- Evidence (`grobench.ts`, `grobench2.ts`): 1crn (327 atoms) × 1000 frames 138–164 ms, 65 MB retained = 200
  B/atom/frame; 13k atoms × 20 frames (`big.gro`, `mkbiggro.mjs`) 98 B/atom/frame. Conformations alone are 12
  B/atom/frame. A 50k-atom, 1000-frame GRO is ~5 GB instead of ~0.6 GB.
- Fix: build one `BasicData` for all frames (concatenate atoms, `pdbx_PDB_model_num` = frame index, as `xyz.ts` does)
  and call `createModels` once; then set `ModelSymmetry` per model from each frame's box. Or build frame 0 and use
  `Model.trajectoryFromModelAndCoordinates` with coordinate frames (as `lammps-trajectory.ts`), but see B7.
- Risk: frames with different topology must still get their own hierarchy (the `sameAsPrevious` check handles it).

### P4: multi-model hierarchy comparison

- Mechanism: `getAtomicHierarchy` runs `findHierarchyOffsets` and `createHierarchyData` (`mapToArray` of
  `type_symbol`, `Table.view` of residues and chains, `columnToArray`, `getResidueSourceIndex` sort) for every model;
  only then `isHierarchyDataEqual` (`:144`) runs `Table.areEqual` on atoms / residues / chains. For text CIF the
  columns are `viewFull` closures, so both sides call `str10` / `int9` and allocate per value. All of it is thrown away
  when the model equals the previous one (NMR, MD).
- Evidence (profile `p-multi50`, `multi.ts`, rerun): multi50.cif (7qpd × 50, 576k rows) 176–186 ms;
  `getAtomicHierarchy` 99 ms, `isHierarchyDataEqual` 65 ms (36 %), `str10` 35 ms self. BCIF (`multib.ts`): compares
  ~10 ms of 132–178 ms.
- Fix: decide `sameAsPrevious` before building: compare the offsets with `previous.atomicHierarchy` segments (Int32),
  then the hierarchy columns row by row on the unsorted `data.atom_site` through
  `col.areValuesEqual(prevSource[i], curSource[i])` (raw-token compare for text, typed compare for BCIF); build only
  when that fails.
- Risk: the compared column set must equal what `createHierarchyData` builds. Chain-operator name substitution
  (`:227`) runs before today's compare; comparing raw source columns is still equivalent since the substitution is a
  function of them. Test multi-model PDB (conversion path), multi-model BCIF, differing atom counts.

### P3: IHM ensembles

- Mechanism: `readIntegrative` calls `createIntegrativeModel` (`:122`) per `ihm_model_list` row: `getCoarse` (keys,
  ranges, index), `getSequence` (all of `entity_poly_seq`), `getAtomicRanges`, `getStructAsymMap` and a fresh
  `_staticPropertyData` (so struct_conn, secondary structure etc. are recomputed per model). Nothing is compared with
  the previous model. Covariances are one `Mat3` (a 9-element JS array) per gaussian
  ([coarse.ts:89](../../packages/model/src/formats/structure/basic/coarse.ts#L89)).
- Evidence (`ihmshare.ts`, `ihmmem2.ts`, rerun): `ihm100.cif` (8zzh × 100; 1764 spheres + 1764 gaussians each)
  225–247 ms, 72 MB retained (~720 KB/model; ~45 KB of it coordinates). Distinct between model 0 and 1: sequence,
  atomicHierarchy, atomicRanges, coarseHierarchy, properties, customProperties, `_staticPropertyData`. For comparison,
  multi50.cif retains 11 MB for 50 × 11.5k atoms.
- Fix: compare `ihm_sphere_obj_site` / `ihm_gaussian_obj_site` (entity_id, asym_id, seq_id_begin, seq_id_end) and the
  atomic part with the previous model through raw-row `areValuesEqual`; when equal reuse coarseHierarchy, sequence,
  atomicHierarchy, atomicRanges, properties, entities and `_staticPropertyData`, and build only the conformation.
  `getSequence` with `entity_poly_seq` depends only on `BasicData`: compute once per file. Store covariances as one
  `Float32Array(9 * count)`.
- Risk: the flat covariance changes `CoarseGaussianConformation`; update graphics (gaussian coarse) and export. Shared
  static props must be topology-only, as for atomic multi-model.

### P8: grid gradients

- Mechanism: `getGradients` computes a `Float32Array(3n)` of central differences through closure `get()` (6 per voxel;
  for `[0,1,2]` layouts `i` innermost, against memory order) and keeps it on `(grid as any)._gradients` for the life of
  the volume. Only consumer: streamlines (`model/src/props/volume/streamlines/basic.ts:167`), which samples along
  seeded lines. The `magnitude` stats it also returns are never read.
- Evidence: 12 B/voxel; 201 MB at 256³, 1.6 GB at 512³ (data itself 537 MB).
- Fix: compute on demand in `makeGetInterpolatedGradient` (central differences of the 8 surrounding cells), nothing
  stored; or keep the precompute on the streamlines property (disposed with it), loop in memory order, drop the unused
  stats.
- Risk: streamlines move slightly if gradients are no longer interpolated from cell values; compare fixtures.
- Addressed (#2017): `makeGetInterpolatedGradient` computes the 8 corners on demand (rounded to float32 as the cache
  stored them: streamlines unchanged); `getGradients` no longer caches. 256³: 201 MB → 0 retained, 532 → 139 ms at the
  default seed density; at the maximum density (30) 3.37 → 3.87 s (+15 %): on-demand reads cost more than the precompute
  when every voxel is sampled many times.

### Low-impact findings

- P7: `getIsosurfaceBoundingSphere` loops over every voxel with `getCoords` per hit; `isobs.ts`: 256³ 35–67 ms (axis
  order 012) and 47–121 ms (210) against 11–17 ms for a plain threshold pass. Called after every CPU mesh, wireframe and
  dot build ([isosurface.ts:205](../../packages/graphics/src/repr/volume/isosurface.ts#L205), `:505`, `dot.ts`); the GPU
  and texture-mesh paths already use `Grid.getBoundingSphere` instead, and CPU marching cubes costs more than the scan.
  Fix: nested loops in memory order (coordinates are loop indices), first / last hit per row, cache by isoValue on the
  grid like the histogram. Addressed (#2017): memory order, first / last hit per row, cached for the last iso value;
  256³ random data 33–122 → 1–6 ms, a blob 34–40 → 8–9 ms (better than the 11–17 ms estimate: rows stop at their first
  hit).
- P5: after bucketing, the identity shortcut needs `indices.length === _rowCount`, so models 2..N and every sorted file
  go through `Table.view` of all schema columns; for BCIF `arrayView` copies each (unused ones included) and
  `Cartn_x/y/z` are copied again into Float32 by `getConformation`. `arrayView` is 13.8 ms of 90 (big-bcif), 23.8 ms of
  137 (multi50.bcif). The 99 MB retained for multi50.bcif (`multib2.ts`, 172 B/atom/model) is the decode of every column
  of the full table (`getNormalizedAtomSite` spreads the lazy `CategoryTable` getters), not the views; per model only id
  / occupancy / B copies stay (~0.15 MB). Fix: `Table.window` when `indices[i] === start + i`; lazy per-column views;
  `getNormalizedAtomSite` replacing only the five substituted columns over the original getters. Addressed (#2017) for
  in-order ranges only: typed columns are windows (no copy), plain-array columns still views (windows of them lost the
  hierarchy compare's array fast path and the time). multi50.bcif 131–134 → 118–122 ms, alloc −35 MB, retained buffers
  66 → 29 MB (the copies were retained, unlike the note above says). Not done: lazy views for sorted models (big-bcif
  63–64 ms unchanged), `getNormalizedAtomSite`.
- P6 (landed, #2007): `${model}|${asym}|${seq}` per atom row (`:83`); same in
  [pdb/entity.ts:235](../../packages/model/src/formats/structure/pdb/entity.ts#L235); `MissingResidues.has/get`
  allocate a key per call. 8.9 ms self on big-bcif. Fix: emit a key only when the triple changes
  (`areValuesEqual(i - 1, i)`), or nested `Map<asym, Map<model, Set<seq>>>` used by `has/get` too.
- P11: `guess.ts`: 2M atoms 100 ms per-atom, 25–27 ms with a call-local `(compId, atomId)` cache. Used by GRO (per
  frame, see P2), PSF, TOP, MOL2, CIF-core. Addressed (#2017): char codes for ASCII names, the regex path for others; 2M
  atoms 83 → 31 ms (no cache needed).
- P9: `new Array(slotCount)` for operatorA/B, distance, flag. Fix: `Int32Array` / `Float32Array` / `Uint8Array` (flags
  fit), or skip constant ones behind `hasOperators`. Addressed (#2017): `Int32Array` operators, `Uint16Array` flags,
  `Float64Array` distances (exact; they are compared with a 0.3 Å tolerance). 1M bonds: −28 MB retained (estimate −36 MB
  assumed Float32 distances); `fromData` 56–61 → 40–42 ms.
- P14: `range(n).sort(comparator)` per model build. Fix: if `residueSourceIndex` is strictly increasing fill 0..n-1.
  Addressed (#2017): 32.4k residues 0.47 → 0.17–0.25 ms (V8's sort was already near-linear on sorted input).
- P13: `Tensor.Space(..., Float64Array)` in `volumeFromCube` / `volumeFromDx`; orbital extraction allocates a new
  Float64 array while `sourceData` keeps the interleaved one. Converting only here adds a copy; the io readers
  (`io/reader/cube`, `io/reader/dx`) must parse to Float32.
- P10: `_p` loops over all entity ids calling `findResidue` (two Map lookups per miss). Resolve the entity from
  `structAsymMap.get(asymId)?.entity_id`, fall back to the loop. Computed once per hierarchy, so small. Addressed
  (#2017): only the entities that have a chain with the asym id, in entity order (same entries, malformed files
  included); 10k rows over 100 entities 37.8 → 4.7 ms.
- P12: `Table.pickRow` (linear) over `ihm_model_group_link` and `ihm_model_group` per model. Build two Maps once.
  Addressed (#2017): `getModelGroupNames(data)` reads both tables once; 5000 models 66 → 0.7 ms, 1000 models 3.3 → 0.5
  ms.

## Bugs

| Id | Title | File | Severity | Verdict |
| --- | --- | --- | --- | --- |
| B2 | Dynamic format properties (`ModelSymmetry`, `IndexPairBonds`) are `undefined` for models 2..N of multi-model files | [property.ts:74](../../packages/model/src/formats/structure/common/property.ts#L74) | high | confirmed |
| B7 | `trajectoryFromModelAndCoordinates` frames drop `IndexPairBonds` set with `Provider.set` (SDF, MOL, MOL2, CIF-core, LAMMPS data) | [model.ts:122](../../packages/model/src/model/structure/model/model.ts#L122) | medium | confirmed |
| B1 | Leak: `getUniqueComponentNames` `memoize1` keeps every `BasicData` (and its source text) | [properties.ts:177](../../packages/model/src/formats/structure/basic/properties.ts#L177) | medium | confirmed |
| B3 | ADP mapping drops the last atom, and every atom when ids are not 1..N | [anisotropic.ts:53](../../packages/model/src/formats/structure/property/anisotropic.ts#L53) | medium | confirmed |
| B4 | `Sequence.ofSequenceRanges` seqIds off by one; caller drops the last coarse element | [sequence.ts:212](../../packages/model/src/model/sequence/sequence.ts#L212) | low | confirmed |
| B5 | PDB CONECT triple bonds written as `'trpl'`, read back as order 1 | [conect.ts:131](../../packages/model/src/formats/structure/pdb/conect.ts#L131) | low | confirmed |
| B6 | `Loci.GranularityOptions` groups plain granularities under "With Symmetry" | [loci.ts:304](../../packages/model/src/model/loci.ts#L304) | low | confirmed |

### B2: dynamic format properties on models 2..N

- Mechanism: `sameAsPrevious` models spread `...previous` ([parser.ts:72](../../packages/model/src/formats/structure/basic/parser.ts#L72)),
  so they share `customProperties` but get a fresh `_dynamicPropertyData`. `FormatPropertyProvider.get` returns early
  on `model.customProperties.has(descriptor)`; once model 0 has computed the value, the shared set has the descriptor
  and model k gets `undefined`. Affects every `asDynamic` provider: `IndexPairBonds` (`molstar_bond_site`) and, since
  `6dfd2b734` (2026-08-08), `ModelSymmetry`.
- Failure (`dynprop.ts`, `dynsym.ts` on `multi2.cif`, 7qpd × 2 with a `_molstar_bond_site` loop): model 0 bonds 2,
  model 1 `undefined`; `ModelSymmetry` model 0 has 1 assembly, model 1 `undefined`. In the plugin, model index ≥ 1 of
  any multi-model mmCIF/PDB gets no assemblies or spacegroup (`root-structure.ts:170/242` fall back to the plain model,
  "has no assembly"), and frames lose explicit bonds. Querying model k first breaks model 0 instead.
- Fix: track presence in the per-model store (`if (name in store) return store[name]`) and skip the
  `customProperties.has` short-circuit for `asDynamic` providers. Shared `customProperties` also means
  `customProperties.assets(desc, …)` for a dynamic custom model property on one model disposes another model's assets
  (`props/common/custom-model-property.ts:76`): plausible, not reproduced.
- Test: two-model `molstar_bond_site` file; `ModelSymmetry` on model 1 before and after model 0.

### B7: coordinate trajectories lose explicitly set bonds (adjacent, `model/structure`)

- Mechanism: `_trajectoryFromModelAndCoordinates` spreads the model with a fresh `_dynamicPropertyData` and re-sets only
  `ModelSymmetry`. Formats that store bonds with `IndexPairBonds.Provider.set` (mol.ts:91, mol2.ts:210,
  lammps-data.ts:155, cif-core.ts:238) have no `formatRegistry` entry, so `get` finds nothing on any frame.
  `trajectoryFromTopologyAndCoordinates` re-sets them; this path does not.
- Failure (`trajbonds.ts`, `long_animation.sdf` + 3 frames): model bonds 11, frames 0–2 `undefined`. SDF / MOL2 + XTC
  / DCD in the plugin (`trajectory-helpers.ts:21`) and LAMMPS data + dump fall back to distance bonding: bond orders
  and non-standard bonds lost on every frame.
- Fix: copy the source model's `IndexPairBonds` (and other dynamic, topology-only props) onto each frame, as the
  topology path does.

### B1: `getUniqueComponentNames` leak

- Mechanism: `memoize1` ([memoize.ts:29](../../packages/core/src/util/memoize.ts#L29)) caches in a module-level `Map`
  keyed by `BasicData`, never evicted. Reached from `getChemicalComponentMap` (`:124`) and
  `getSaccharideComponentMap` (`:169`) whenever `chem_comp` is empty. The key holds `atom_site`, whose CIF field closures
  hold the tokens and the source string.
- Failure (`leak2.ts`, `multi2.cif`, distinct string per load): heap +2.7 MB per load, never released; 00-profile:
  +4.1 MB per load of 7qpd.fw2.cif. Unbounded in long sessions for any mmCIF without `chem_comp`.
- Fix: compute `getUniqueComponentNames(data)` once in `getCommonProperties` when `chem_comp` is empty and pass it to
  both; or a `WeakMap`.

### B3: anisotropic ADP mapping

- Mechanism: `atomIdToElement = new Int32Array(atomId.rowCount)` indexed by the id value: 1-based id N writes past the
  end; ids > rowCount (models 2..N, PDB serials with TER gaps, subsets) are dropped; the lookup then reads `undefined`,
  and `undefined !== -1` writes to `elementToAnsiotrop[undefined]` (no-op).
- Failure (`aniso.ts`): ids 1..3 with ADPs for all → `[0, 1, -1]`; ids 101..103 → `[-1, -1, -1]`. Ellipsoids missing.
  Used by mmCIF and PDB ANISOU (through `pdbToMmCif`).
- Fix: a `Map<id, element>` (or an `Int32Array` sized `max(id) + 1` when reasonable); test `ei !== undefined`.

### B4–B6

- B4: `seqId.value(row) = row + minSeqId + 1` but `index(seqId) = seqId - minSeqId` (`seqranges.ts`: ranges [1,11],
  [10,20] → seqId(0) = 2, last = 21). Also
  [properties/sequence.ts:131](../../packages/model/src/model/structure/model/properties/sequence.ts#L131) passes
  `offsets[cI + 1] - 1` as the exclusive `Column.window` end, dropping the last coarse element per entity. Only coarse
  models without `entity_poly_seq` (`StructureSequence.fromHierarchy`); 8zzh has `entity_poly_seq` and is unaffected
  (`coarseseq.ts`). Consumers: sequence viewer, `AlignSequences.compute`, atomic ranges. Fix: `row + minSeqId` and
  `offsets[cI + 1]`.
- B5: `['sing', 'doub', 'trpl', 'quad']`; struct_conn reading expects `'trip'`
  ([struct_conn.ts:167](../../packages/model/src/formats/structure/property/bonds/struct_conn.ts#L167)), so
  `CONECT 1 2 2 2` becomes order 1 (ligand table fallback). Fix: `'trip'`.
- B6: `k.indexOf('Instances') ? …` is truthy for -1; residue / chain / entity / model / operator get "With Symmetry"
  (`gran.ts`). Fix: `k.includes('Instances')`.

## Checked and dropped

Nothing refuted. Corrections to the audit:

- P5: the memory claim ("several MB less retained per BCIF model") is wrong; the retained bulk is the full-table BCIF
  decode, needed for the used columns. Downgraded to low.
- P7: the GPU isosurface paths already avoid the scan; only CPU mesh, wireframe and dots pay it. Downgraded to low.
- P2: 200 B/atom/frame holds for small frames; at 13k atoms it is ~100 B/atom/frame (the 10 GB extrapolation is
  ~5 GB).
- P1: flat Int32 gives 2.3–5× in time (5× was the best case).
- P4 risk note: chain-operator substitution runs before the compare, not after.
- Memory measured in a loop with a "previous run" baseline is skewed by `memoizeLatest` caches and by B1 (`run 1`
  numbers ≈ 0 or negative); use one fresh process per number (`grobench2.ts`, `ihmmem2.ts`, `multib2.ts`).

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/model-formats/`,
bundled with `node ../build.mjs x.ts x.mjs`, run with `node --expose-gc x.mjs <file>`:

| Script | Finding | Input |
| --- | --- | --- |
| `align.ts`, `align2.ts` | P1 | synthetic sequences |
| `pdbbench.ts` | P1 | `seqres.pdb` (`mkseqres.py`), `noseqres.pdb` |
| `grobench.ts`, `grobench2.ts`, `mkbiggro.mjs` | P2 | `multi1000.gro`, `multi200.gro`, `big.gro` |
| `ihmshare.ts`, `ihmmem2.ts` | P3 | `ihm100.cif` (`mkihm.py`) |
| `multi.ts`, `multib.ts`, `multib2.ts`, `../analyze2.mjs` (`MAP=multi.mjs.map`) | P4, P5 | `multi50.cif`, `multi50.bcif` (`mkmulti.py`, `tobcif.ts`), profiles `p-multi50*` |
| `isobs.ts` | P7 | synthetic 256³ grids |
| `guess.ts` | P11 | synthetic 2M atoms |
| `dynprop.ts`, `dynsym.ts` | B2 | `multi2.cif` |
| `trajbonds.ts` | B7 | `data/examples/long_animation.sdf` |
| `leak2.ts`, `leak1.ts` | B1 | `multi2.cif` |
| `aniso.ts`, `seqranges.ts`, `coarseseq.ts`, `gran.ts` | B3, B4, B6 | synthetic / 8zzh |
