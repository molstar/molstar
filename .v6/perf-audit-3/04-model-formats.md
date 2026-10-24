# Slice 04: model formats and model (v3)

Scope: [`packages/model/src/formats`](../../packages/model/src/formats) (structure, volume, shape and topology
converters, format properties, writers), [`packages/model/src/model`](../../packages/model/src/model) without
`model/structure` (sequence, alignment, volume grid, loci, custom properties), and the kinemage reader
([`packages/extensions/kinemage/src/reader`](../../packages/extensions/kinemage/src/reader)). `model/periodic.ts` and
`Model.withAtomicConformation` (#2037) are reviewed in [05](05-model-structure.md); this slice covers #2037's reader
side (GRO, LAMMPS, DCD, NAMD XSC). HEAD `16277d822`, v2 base `0944823d4`, Node 26.12, Apple M4. Read-only.

Measurement conditions: machine overloaded throughout (load 15–38, free memory ~60 MB, ~13 other auditors). Wall
times swung 2–5× between identical runs (flat-192 `getMissingResidues` alone: 170–1000 ms for a ~20 ms loop), so the
claims below rest on **allocation** (`total_allocated_bytes`) and **retained** bytes (heap + array buffers after GC),
which reproduce to ±2 %; wall times are ranges from base/HEAD runs interleaved back-to-back. A/B base:
`git archive 0944823d4` in `pa3/04/base`, bundled by the same esbuild script.

PRs since v2 touching the slice: #2027 (kinemage moved to `packages/extensions`), #2032 (B7, 12/B3), #2033 (B3, B5,
B6, B11), #2036 (P5, P16 part: lazy sorted columns, Float32 gather, dictionary string columns, bucket keys on
dictionary indices), #2037 (periodic: GRO triclinic box, LAMMPS boxes/tilt/image flags, DCD step, NAMD XSC, new
code), #2039 (text CIF from bytes, `Sniffable` in `mmcif-sniff`), #2044 (B4, B12, B13), #2045 (B14), #2047 (P6, P17
part, P18, 12/B1, 12/L2 part). #2016 (bond providers) touched `struct_conn` order flags.

## Hot paths now

| Path | Cost (fixture) | Where it goes |
| --- | --- | --- |
| single-model BCIF, 1.12M atoms (flat-192), model stage | 0.63–0.87 s wall loaded (base 1.0–1.33 s same window); **112–179 MB alloc (base 487–550)**; retained **124 B/atom = 139 MB (base 320 = 359 MB)** | profile (loaded, 1363 ms sampled): BCIF decode moved here lazily (`fusedNumbers` 14 %, `decodeStringDictionary` 6 % self: io); `addMissingResiduesFromEntityPolySeq` 10 % self (P15); `Column.view` of sorted columns 17 % total; `getAtomicHierarchy` 37 % total; `_makeBuckets` 7 % |
| flat-32 (186k atoms) model stage | 155–424 ms (base 168–399, noise-bound); **25–40 MB alloc (base 84–98)**; retained 125 B/atom (base 323) | as flat-192 |
| multi-model text mmCIF, 50 × 11.5k (multi50) | 427–484 ms (base 497–579), **42–47 MB alloc (base 153–159)** | P4: hierarchy built per model then compared; `str` per row 21 %, byte-string cache `read` 20 %, `areValuesEqual` 8 % self |
| multi-frame GRO, 1000 × 327 atoms | 184 B/atom/frame retained (base 178) | P2: a full model per frame |
| LAMMPS dump with image flags, 10k atoms × 100 frames | 24.3 B/atom/frame retained (12.7 without images) | P20: unwrapped copies made eagerly for every frame |
| NAMD XSC + DCD, 30k frames × 30k entries | 665 ms (loaded) | P21: per-frame linear scan of the entries |
| Kinemage 500k dots + 500k vectors | 208 B/point (base 316) | P17: per-point position arrays and labels |

Retained breakdown of a sorted BCIF model (flat-32 heap snapshot, 186k atoms, `snap.ts` + `snapan.mjs`): source
decoded columns ~60 B/atom (Int32 20, Float64 occupancy/B 16, Float32 Cartn cache 12, dictionary indices 10, masks 2),
model copies ~65 B/atom (Int32 atomId/sourceIndex/segments/charge ~24, Float64 occupancy/B views 16, Float32 x/y/z 12,
`type_symbol` JS array 8, dictionary indices ~5), encoded bytes 16 B/atom (io). See P16.

## v2 → now

28 v2 items (14 perf, 14 bugs): **fixed 15, partly fixed 4, still open 9, regressed 0, dropped 0**. Every probe
rerun at HEAD (`v1re/*`, `srcidx*`, `conect`, `stats`, `missing`, `kin`, `model`, `leak`).

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P5 | Sorted models copy every `atom_site` column; normalization decodes every column | fixed | #2036 | lazy columns ([sort.ts:109](../../packages/model/src/formats/structure/basic/sort.ts#L109)), BCIF fields decode on first read, Cartn gathered to Float32; flat-192 alloc 487–550 → 112–179 MB |
| P16 | A BCIF model retains ~320 B/atom | partly fixed | #2036 | 320 → 124 B/atom (flat-192 359 → 139 MB); remainder below |
| P4 | Multi-model mmCIF builds each hierarchy before comparing | still open | | [atomic.ts:241](../../packages/model/src/formats/structure/basic/atomic.ts#L241) unchanged; multi50 427–484 ms, alloc 42–47 MB (allocation −70 % via #2036/#2039, compare still per row through `str`) |
| P2 | Multi-frame GRO / MOL2 build a full model per frame | still open | | [gro.ts:139](../../packages/model/src/formats/structure/gro.ts#L139); 184 B/atom/frame (v2 200; rerun base 178) |
| P3 | IHM ensembles rebuild everything per model | still open | | [parser.ts:280](../../packages/model/src/formats/structure/basic/parser.ts#L280) and `coarse.ts` unchanged since v2 (not rerun: `ihm100.cif` gone with the v1 scratchpad) |
| P15 | `missingResidues` computed eagerly on every load | still open | | [parser.ts:49](../../packages/model/src/formats/structure/basic/parser.ts#L49); now 10 % self of the flat-192 model stage (the rest shrank) |
| P17 | Kinemage: tag array per point, JS arrays for positions | partly fixed | #2047 | rerun 316 → 208 B/point; per-point `position` array and label open ([kinparser.ts:118](../../packages/extensions/kinemage/src/reader/kinparser.ts#L118)) |
| P19 | `molstar_bond_site` bonds rebuilt per model and per coordinate frame | partly fixed | #2032 | frames now copy the topology model's `IndexPairBonds` ([model.ts:117](../../packages/model/src/model/structure/model/model.ts#L117)); per model of a multi-model file still `toDatabase` + all rows ([mmcif.ts:95](../../packages/model/src/formats/structure/mmcif.ts#L95)) |
| P13 | Cube and DX volumes as Float64 | still open | | [cube.ts:21](../../packages/model/src/formats/volume/cube.ts#L21), [dx.ts:18](../../packages/model/src/formats/volume/dx.ts#L18) |
| 12/T1 | Trilinear interpolation: 8 `get()` closure calls per sample | still open | | [grid.ts:214](../../packages/model/src/model/volume/grid.ts#L214); #2032 added 6 modulo ops per periodic sample |
| P6 | PDB SEQRES scan: key per atom | fixed | #2047 | key on change ([entity.ts:240](../../packages/model/src/formats/structure/pdb/entity.ts#L240)) |
| P18 | Volume stats: four passes | fixed | #2047 | core `arrayStats`; `stats.ts` 256³: 116–165 → 39–91 ms (loaded); `stats.test.ts` counts reads |
| 12/B1 | Segmentation bounds: `getCoords` per voxel | fixed | #2047 | memory-order loops ([segmentation.ts:100](../../packages/model/src/formats/volume/segmentation.ts#L100)) |
| 12/L2 | `Cell.getBoundingSphere` and CCP4 writer: accessor per cell | partly fixed | #2047 | writer by strides (row `set`); [volume.ts:351](../../packages/model/src/model/volume/volume.ts#L351) still two passes, `getCoords` + 2 transforms per cell per instance |
| B8 | `molstar_bond_site` bonds in element order, read as source rows | still open | | `srcidx2.ts`: 0 of 2 bonds; `srcidx.ts`: model 2 IndexPairBonds 3 edges, unit bonds 0 ([mmcif.ts:106](../../packages/model/src/formats/structure/mmcif.ts#L106)) |
| B10 | Volume sigma is RMS, not standard deviation | still open | | #2047 kept `rms` in `arrayStats` ([ccp4.ts:132](../../packages/model/src/formats/volume/ccp4.ts#L132), cube, dx, dsn6); `stats.ts`: sigma 103.26 for sd 7.35 |
| B1 | `getUniqueComponentNames` `memoize1` keeps every `BasicData` | still open | | [properties.ts:189](../../packages/model/src/formats/structure/basic/properties.ts#L189); `leak.ts` 7qpd (no `chem_comp`): +3.0 MB per load, 8 loads +25.3 MB (base +2.5 MB/load; HEAD keeps the UTF-8 bytes of a string input) |
| B7 | Coordinate frames lose `IndexPairBonds` set with `Provider.set` | fixed | #2032 | `trajbonds.ts`: model 11, frames 11 |
| B3 | ADP mapping drops the last atom / ids not 1..N; shared models | fixed | #2033 | `aniso.ts`: `[0,1,2]` for both; provider `asDynamic` |
| 12/B3 | Periodic trilinear interpolation reads unwrapped indices | fixed | #2032 | floor and wrap both neighbours ([grid.ts:251](../../packages/model/src/model/volume/grid.ts#L251)) |
| B9 | Models sharing the previous hierarchy get its `atomSourceIndex` | still open | | `srcidx.ts`: model 2 source index `[0,1,2]` (its rows are 3..5); [atomic.ts:243](../../packages/model/src/formats/structure/basic/atomic.ts#L243) |
| B11 | One-sided CONECT dropped | fixed | #2033 | `conect.ts`: 1 entry for each of the three cases |
| B4 | `Sequence.ofSequenceRanges` off by one | fixed | #2044 | `seqranges.ts`: seqId(0) 1, last 20 |
| B5 | CONECT triple bonds as `'trpl'` | fixed | #2033 | |
| B6 | Granularity options under 'With Symmetry' | fixed | #2033 | `gran.ts`: only `*Instances` grouped |
| B12 | MOL2 bonds assume ids 1..N | fixed | #2044 | id → row map, unknown ids dropped |
| B13 | Segmentation throws on a label without a set | fixed | #2044 | `if (!b) continue` in both loops |
| B14 | Kinemage re-opened `@group` keeps previous masters | fixed | #2045 | `groupMastersByName` |

New costs or bugs from fixes since v2: **P20** and **P21** (new code in #2037). No regression: the model stage
allocates 65–77 % less and retains 61 % less than base on flat-32/192; multi50 −70 % alloc, no slower. The GRO parse
result retains 21 MB vs 6 MB in base only because the bench passes a string, which the byte tokenizer (#2040) encodes
and keeps; the app passes bytes (io, not this slice). One pre-existing bug found while reading the #2036 hunk (B15).

## Open performance findings (ranked)

| Id | Status | Title | File | Kind | Estimate | Effort | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P16 | partly | BCIF model still retains 124 B/atom: occupancy/B as Float64 twice, Cartn twice on sorted models, `type_symbol` as a JS array | [atomic.ts:145](../../packages/model/src/formats/structure/basic/atomic.ts#L145) | memory | −40 B/atom estimated from the snapshot: −45 MB at 1.1M, −200 MB at 5M | S–M | no |
| P4 | still open | Multi-model mmCIF builds each hierarchy before comparing | [atomic.ts:241](../../packages/model/src/formats/structure/basic/atomic.ts#L241) | time | ~450 → ~250 ms per 576k rows (estimated: compare + build are ~45 % of the profile) | M | no |
| P2 | still open | Multi-frame GRO / MOL2: full model per frame | [gro.ts:139](../../packages/model/src/formats/structure/gro.ts#L139) | memory | 184 → ~15–25 B/atom/frame (v2 estimate) | S | no |
| P20 | new | LAMMPS: unwrapped coordinates made eagerly for every frame | [lammps-trajectory.ts:111](../../packages/model/src/formats/structure/lammps-trajectory.ts#L111) | memory | measured: 24.3 → 12.7 B/atom/frame; 100k atoms × 1000 frames −1.2 GB | S | no |
| P3 | still open | IHM ensembles rebuild per model | [parser.ts:280](../../packages/model/src/formats/structure/basic/parser.ts#L280) | time+memory | ~720 → ~60–100 KB/model (v2) | M | no |
| P15 | still open | `missingResidues` eager | [parser.ts:49](../../packages/model/src/formats/structure/basic/parser.ts#L49) | time | 10 % of the flat-192 model stage; −18–31 ms at 1.1M (v2, idle) | S | no |
| P21 | new | NAMD XSC: each frame scans the entries from the start | [namd-xsc.ts:36](../../packages/model/src/formats/structure/namd-xsc.ts#L36) | time | measured (loaded): 1k 12 ms, 10k 109 ms, 30k 665 ms; O(F·E) → O(F + E) | S | no |
| P17 | partly | Kinemage per-point position arrays and labels | [kinparser.ts:118](../../packages/extensions/kinemage/src/reader/kinparser.ts#L118) | memory | 208 → ~60–80 B/point (v2) | S | no |
| P19 | partly | `molstar_bond_site` per model of a multi-model file | [mmcif.ts:95](../../packages/model/src/formats/structure/mmcif.ts#L95) | time | O(M·(A + M·B)) → O(A + M·B) (by code) | S | no |
| P13 | still open | Cube and DX as Float64 | [cube.ts:21](../../packages/model/src/formats/volume/cube.ts#L21) | memory | −50 % volume memory | S | no |
| 12/T1 | still open | Trilinear interpolation through a closure, 8 calls/sample | [grid.ts:214](../../packages/model/src/model/volume/grid.ts#L214) | time | ~1.4× per sample (v2, plausible) | S | no |
| 12/L2 | partly | `Volume.Cell.getBoundingSphere`: two passes, accessor per cell | [volume.ts:351](../../packages/model/src/model/volume/volume.ts#L351) | time | ~2× per call | S | no |

### P16: what a BCIF model still retains

- Evidence (`snap.ts` heap snapshot on flat-32, grouped by retainer in `snapan.mjs`; per atom):
  - occupancy and B: the source decodes Float64 (16 B), and the sorted model keeps `Column.view` Float64 copies (16 B)
    through `floatRows` ([sort.ts:134](../../packages/model/src/formats/structure/basic/sort.ts#L134)) because
    `getConformation` stores the columns, read later per atom by themes;
  - Cartn: the field's Float32 decode is cached (`binary/field.ts` `float32`, 12 B) and the sorted gather copies it
    (12 B). In-order models share the cache (window), so this half is sorted-only (9rxg, 21nr, Mol*-written
    assemblies, the fixtures);
  - `atoms.type_symbol`: `Column.mapToArray(…, ElementSymbol)` into a JS array, 8 B/atom (4 in Chrome)
    ([atomic.ts:66](../../packages/model/src/formats/structure/basic/atomic.ts#L66)), while the source is a dictionary.
- Fix: (1) `getConformation` stores `Column.ofFloatArray(occupancy.toArray({ array: Float32Array }))` and the same for
  B: a quantized BCIF field decodes straight to Float32 and the Float64 decode never happens (−24 B/atom); (2) the
  sorted gather decodes without filling the field's Float32 cache (an uncached `toFloatArray`, io) (−12 B/atom); (3)
  `type_symbol` as a dictionary: map the dictionary values through `ElementSymbol` once, keep the indices (−7 B/atom).
- Risk: (1) changes the B/occupancy column type seen by exports (written with 2 decimals; Float32 holds them); keep
  `isDefined`/value kinds.

### P20: LAMMPS unwrapped coordinates

- `lammpsFrameCoordinates` builds `ux, uy, uz` for every frame with image flags at load; `model.ts:143` already makes
  the per-frame `UnwrappedConformation` lazily, but over arrays that exist for all frames. `lmp.ts`, 10k atoms × 100
  frames: retained 12.7 → 24.3 B/atom/frame, alloc 23.5 → 42.9 MB, conversion time +33 % (loaded).
- Fix: `Frame.unwrapped` as a thunk (or the image columns plus the box) and compute it in the per-frame getter;
  the parsed file keeps `frame.image` anyway.

### P21: XSC entry per frame

- `entryFor` loops `for (const e of entries)` from the start for each frame. `xsc.ts` (frames = entries, steps
  every 500): 1k 12 ms, 10k 109 ms, 30k 665 ms (loaded), quadratic; a 100k-frame NAMD run ~7 s at load time.
- Fix: frames and entries are both in step order: one advancing pointer (or binary search on `step`).

### P4, P2, P3, P15, P17, P19, P13, 12/T1, 12/L2

Unchanged code and fixes as in v2 ([../perf-audit-2/04-model-formats.md](../perf-audit-2/04-model-formats.md)).
P4's profile now: the comparison (`areColumnsEqual` → `str`, the byte-string cache `read`) is ~45 % of the stage, so
"compare source rows with the previous model's before building" still halves it.

## Bugs

| Id | Severity | Status | Title | File | Failure scenario |
| --- | --- | --- | --- | --- | --- |
| B8 | medium | still open | `molstar_bond_site` bonds built in element order, read as source rows | [mmcif.ts:106](../../packages/model/src/formats/structure/mmcif.ts#L106) | reload a sorted Mol*-written mmCIF: 0 of 2 bonds; model 2 with its own topology: 3 IndexPairBonds edges, 0 unit bonds |
| B10 | medium | still open | Volume sigma is RMS, not SD (CCP4 without ARMS, DX, Cube, DSN6) | [ccp4.ts:132](../../packages/model/src/formats/volume/ccp4.ts#L132) | mean 103, sd 7.35, max 117: sigma 103.26, 1σ above the max, empty default isosurface |
| B1 | medium | still open | `memoize1` in `getUniqueComponentNames` keeps every `BasicData` | [properties.ts:189](../../packages/model/src/formats/structure/basic/properties.ts#L189) | each mmCIF load without `chem_comp` keeps ~3 MB (7qpd) forever: 8 loads +25.3 MB |
| B9 | low | still open | Shared hierarchy reuses model 1's `atomSourceIndex` | [atomic.ts:243](../../packages/model/src/formats/structure/basic/atomic.ts#L243) | NMR export of model 2 writes model 1's `pdbx_label_index`/SIFTS rows |
| B15 | low | new (pre-existing) | `sortAtomSite` tests `label_seq_id` presence at a bucket position, not at its row | [sort.ts:37](../../packages/model/src/formats/structure/basic/sort.ts#L37) | chain A protein, a water chain, chain B with residues 2 then 1: B is not sorted (`ALA1 GLY2 ALA1 HOH0`; without the water `ALA1 ALA1 GLY2`) |

### B15: seq-id sort skipped by position

- `aI = chainBuckets[cI]` is a position in `indices`; `label_seq_id.valueKind(aI)` reads row `aI`. After entity
  bucketing moved rows, or for any model after the first (`indices` start at the model's first row, positions at 0),
  that row belongs to another chain. A polymer chain whose position maps to a HETATM row is left unsorted; a model k > 1
  tests model 1's rows. `seqsort.ts`: chain B residues stay `GLY2 ALA1`; residue order then disagrees with
  `label_seq_id` (polymer gaps, sequence mapping). Needs out-of-order residues in the file, so rare.
- Fix: `label_seq_id.valueKind(indices[aI])`. Pin with the `seqsort.ts` case in `sort.test.ts`.

### B8, B10, B1, B9

Mechanisms and fixes as in v2. B10: #2047 replaced the four passes with `arrayStats` but kept its `rms`; return
`sd = sqrt(sumSq/n − mean²)` (or Welford) and use it as `sigma` in the four readers. B1: a `WeakMap`, or compute the
names in `getCommonProperties`.

## Already good, re-checked

- #2036 sorted path: `floatRows` gives Float32 for quantized and text sources and keeps the source type for a raw
  Float64 column, as the in-order path; `bucketKey` canonicalizes the duplicate `''` of BinaryCIF dictionaries and
  skips columns with a `transform`; `lazyTable` setters keep `getNormalizedAtomSite`'s substitutions.
- #2037 readers: GRO box vector order (v1x v2y v3z v1y v1z v2x v2z v3x v3y) maps to a, b, c correctly; LAMMPS
  tilt/scaled/`xu` handling and the periodicity guard; DCD `step`; XSC empty files are a parse error, frames before the
  first entry take the first.
- #2033 CONECT: per-side counts give the right order for one- and two-sided listings; `trip`.
- #2047 segmentation strided loop, CCP4 writer row copy (typed `data` only), `arrayStats`; #2032 periodic
  interpolation; #2044 MOL2 id map; anisotropic `Map`.
- Retention: a loaded model and its parse are fully collected once dropped (WeakRefs, `leak2.ts`), except B1.

## Not audited

- `model/periodic.ts`, `Model.withAtomicConformation` and the per-frame model cost: slice 05.
- BCIF decode (`fusedNumbers`, `decodeStringDictionary`), the byte tokenizers and parse retention: io.
- `areBackboneConnected` allocates 4 MB (20 % of the sampled JS allocation) in the flat-32 model stage
  (`atomic-ranges.ts:18`, model/structure): left to 05.
- IHM (P3) not rerun: its input is gone; code unchanged.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/04/`:
`node build.mjs <x.ts> <x.mjs> [base]` (esbuild, `molstar-src`; `base` bundles `base/` = `0944823d4`), run with
`node --expose-gc`. `model.ts` (model stage time/alloc/retained), `leak.ts`, `leak2.ts` (B1, retention), `snap.ts` +
`snapan.mjs` (P16 breakdown), `lmp.ts` (P20), `xsc.ts` (P21), `seqsort.ts` (B15), `kin.ts`/`kinb.ts` (P17),
`srcidx*.ts`, `conect.ts`, `stats.ts`, `missing.ts`, `v1re/*` (P2, P4, B3, B4, B6, B7); inputs `in/multi50.cif`,
`in/multi1000.gro` from `gen.py`; profiles `prof192/`, `profmulti*/`, `alloc32/`.
