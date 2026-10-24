# Performance audit v2: bugs

Every open bug found or re-verified by the v2 slice audits ([README.md](README.md)), confirmed or plausible after the
verify pass, ranked by severity, then by reach (reachable in the apps before latent API defects). 30 distinct entries
(01/B1 = 04/B1 in one row; slice 12's re-listing of 03/B7 merged into the owning slice): **0 high, 10 medium, 20
low**. 14 new (0 high, 7 medium, 7 low), 16 v1 bugs still open (3 medium, 13 low). 3 plausible (05/B3, 05/B13,
07/B9). Ids are `slice/id`; the slice report's Bugs section has the probe and the
fix.

## High

None open: 08/B6 is fixed (below).

## Medium

| Id           | Status vs v1 | Title                                                                                          | File                                                                                                                                                   | Failure scenario                                                                                                                                                       | Verdict   | Report                                                                         |
| ------------ | ------------ | ---------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ------------------------------------------------------------------------------ |
| 05/B14       | new          | `within` without min-radius adds the query atom's vdW radius, ignores `atom-radius`            | [filters.ts:189](../../packages/model/src/model/structure/query/queries/filters.ts#L189)                                                               | PyMOL `around`/`within`/`near_to`/`beyond` over-select by one vdW radius (7qpd, 3 Å: 126 atoms vs 42 by centre distance); `gap` misses the target radius (126 vs 214) | confirmed | [05](05-model-structure.md#bugs)                                               |
| 08/B8        | new          | Smoothed `volumeInstance` colors and layers not recomputed on a transform change               | [visual.ts:744](../../packages/graphics/src/repr/visual.ts#L744)                                                                                       | Overpaint/transparency/emissive/substance on a surface, then move it (repr transform, unitTransforms, assembly unwind): the layer grid is sampled at the old place      | confirmed | [08](08-graphics-geo.md#b8-bug-smoothed-volumeinstance-grids-go-stale-on-a-transform-change) |
| 04/B8        | new          | `molstar_bond_site` bonds built in element order but read as source rows                       | [mmcif.ts:106](../../packages/model/src/formats/structure/mmcif.ts#L106)                                                                               | Reload a sorted Mol\*-written mmCIF (assembly with ligands): 0 of 2 inter-unit bonds; models 2..N with their own topology get no explicit bonds                      | confirmed | [04](04-model-formats.md#b8-molstar_bond_site-index-space)                     |
| 04/B10       | new          | Volume sigma computed as RMS, not standard deviation (CCP4 without ARMS, DX, Cube, DSN6)       | [ccp4.ts:122](../../packages/model/src/formats/volume/ccp4.ts#L122)                                                                                    | Map with mean 103, sd 7.35, max 117: sigma 103.26, so a 1σ level (206) lies above the map maximum and nothing shows                                                  | confirmed | [04](04-model-formats.md#b10-sigma-as-rms)                                     |
| 01/B1, 04/B1 | still open   | `memoize1` never evicts; `getUniqueComponentNames` keyed by `BasicData` pins each load          | [properties.ts:189](../../packages/model/src/formats/structure/basic/properties.ts#L189), [memoize.ts:29](../../packages/core/src/util/memoize.ts#L29) | Every mmCIF load without `chem_comp` keeps its source text: +2.7–4.1 MB per load, forever; plugin/ui `parameters.tsx:962/964` add one entry per color array          | confirmed | [01](01-core-data-util.md#bugs), [04](04-model-formats.md#b1-b7-b3-12b3)        |
| 10/B4        | new          | While any animation (camera spin included) or a snapshot plays, outside edits never make checkpoints unsafe | [checkpoint.ts:205](../../packages/prism/viewer/src/sources/checkpoint.ts#L205)                                                           | A script commits during camera spin, the user presses Undo: the outside edit is silently lost                                                                         | confirmed | [10](10-prism-viewer.md#bugs)                                                  |
| 03/B6        | still open   | `doError` clears objects and caches without transformer `dispose`                              | [state.ts:1257](../../packages/core/src/state/state.ts#L1257)                                                                                          | A structure update returns Null or throws: theme and custom-property refcounts, component descriptors and Download assets below it leak, per re-apply                  | confirmed | [03](03-core-state-task.md#bugs)                                               |
| 03/B11       | new          | Flags carried into a builder commit override a tree passed in after them (from the 03/B1 fix)  | [state.ts:467](../../packages/core/src/state/state.ts#L467)                                                                                            | Builder outstanding across a queued `setSnapshot` or undo: cell hidden-flag undefined, tree true; the saved session and undo history hide what is shown               | confirmed | [03](03-core-state-task.md#bugs)                                               |
| 06/B2        | still open   | Interaction features cached on shared Units but computed from the local substructure's bonds   | [interactions.ts:406](../../packages/model/src/props/computed/interactions/interactions.ts#L406)                                                       | Features depend on which component computed first (7qpd units 1 and 3 differ by one feature)                                                                          | confirmed | [06](06-model-props-script.md#bugs)                                            |

## Low

Reachable today first, then latent (no in-repo trigger; the latent one-liners were fixed in #2046–#2048).

| Id       | Status vs v1 | Title                                                                                        | File                                                                                                                                         | Failure scenario                                                                                                                                  | Verdict   | Report                                                                                |
| -------- | ------------ | -------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ------------------------------------------------------------------------------------- |
| 05/B7    | still open   | `Bundle.toLoci` / `toStructure` create an Error but never throw it                           | [bundle.ts:156](../../packages/model/src/model/structure/structure/element/bundle.ts#L156)                                                   | A bundle outliving its structure silently selects wrong atoms                                                                                     | confirmed | [05](05-model-structure.md#bugs)                                                      |
| 05/B13   | new          | Static component remap trusts the weak (size, min, max, mid) element hash                    | [structure.ts:1111](../../packages/model/src/model/structure/structure/structure.ts#L1111)                                                   | A static component below a source whose elements change between frames ([1,3,5,9] vs [1,4,5,9]) stays stale                                       | plausible | [05](05-model-structure.md#bugs)                                                      |
| 05/B3    | still open   | `Bundle.fromLoci` / `Loci.toExpression` group by a 32-bit hash without equality              | [bundle.ts:119](../../packages/model/src/model/structure/structure/element/bundle.ts#L119)                                                   | Rare wrong atoms in stored selections (P ≈ n²/2³³)                                                                                                | plausible | [05](05-model-structure.md#bugs)                                                      |
| 04/B9    | new          | Models sharing the previous hierarchy also get its `atomSourceIndex` (model 1's rows)        | [atomic.ts:244](../../packages/model/src/formats/structure/basic/atomic.ts#L244)                                                             | NMR/MD multi-model export: model 2's `pdbx_label_index` and SIFTS values come from model 1                                                        | confirmed | [04](04-model-formats.md#low-bugs)                                                    |
| 06/N2    | still open   | Refiners write per-member decisions into a group's shared intra-contact flags                | [refine.ts:190](../../packages/model/src/props/computed/interactions/refine.ts#L190)                                                         | In a symmetry range, a weak H-bond hidden by one copy's crystal contact disappears in every copy                                                  | confirmed | [06](06-model-props-script.md#bugs)                                                   |
| 06/B10   | still open   | Streamline tracer's `prev` point is module state, never reset per line                       | [basic.ts:42](../../packages/model/src/props/volume/streamlines/basic.ts#L42)                                                                | Occasional missing first point of a streamline                                                                                                    | confirmed | [06](06-model-props-script.md#bugs)                                                   |
| 07/B4    | still open   | Segment voxel loci use bbox-local group ids as global cell indices                           | [segment.ts:111](../../packages/graphics/src/repr/volume/segment.ts#L111)                                                                    | 'voxel' granularity on segmentations: the wrong voxel picked and highlighted                                                                      | confirmed | [07](07-graphics-repr.md#low-1)                                                       |
| 07/B8    | still open   | Highlight removal on reused visuals mutates `tMarker` without a cell update                  | [units-representation.ts:169](../../packages/graphics/src/repr/structure/units-representation.ts#L169)                                       | Structure change while hovering: stale highlight until the next mark                                                                              | confirmed | [07](07-graphics-repr.md#low-1)                                                       |
| 07/B9    | new          | `Visual.mark` treats overlapping intervals as full coverage                                  | [visual.ts:215](../../packages/graphics/src/repr/visual.ts#L215)                                                                             | Superposition panel atom entries (two atoms in one residue) on a polymer-granular visual: whole visual highlighted, or a toggle cancelled           | plausible | [07](07-graphics-repr.md#b9-overlapping-intervals-taken-as-full-coverage-new-plausible) |
| 08/B1    | still open   | `Mesh.fillEdges` never adds a triangle (`< 2` filter, then reads `bnd[1]`)                   | [mesh.ts:566](../../packages/graphics/src/geo/geometry/mesh/mesh.ts#L566)                                                                    | `includeParent` surfaces keep their edge holes and notches                                                                                        | confirmed | [08](08-graphics-geo.md#bugs)                                                         |
| 12/B4    | new          | Bench stages swallow errors into `info`, so a failing stage reads as fast                    | [stage.ts:99](../../tests/src/pipeline/stage.ts#L99)                                                                                         | A broken stage shows as a speedup in any comparison (also `bench.ts:209`); blocks 12/H10                                                          | confirmed | [12](12-app-level.md#bugs)                                                            |
| 12/B1    | new          | Pipeline bench crashes on coarse-only entries (`selections()`)                               | [bench.ts:324](../../tests/src/pipeline/bench.ts#L324)                                                                                       | `run --fixture data/examples/ihm-coarse-8zzh.cif`: TypeError reading 'elements' after `theme-eval`; no `--skip` avoids it                         | confirmed | [12](12-app-level.md#bugs)                                                            |
| 01/B10   | still open   | `ungzip` trusts the last member's ISIZE                                                      | [zip.ts:241](../../packages/core/src/util/zip/zip.ts#L241)                                                                                   | `gzip -c a b > ab.gz` (valid multi-member) fails to open; > 4 GiB output likewise                                                                 | confirmed | [01](01-core-data-util.md#bugs)                                                       |
| 01/B14   | new          | number-parser: `parseInt` wraps above 2³¹−1; `parseFloat` not correctly rounded              | [number-parser.ts:39](../../packages/core/src/util/number-parser.ts#L39)                                                                     | '3000000000' → −1294967296 silently; 149 of 200k 3-decimal values differ from `Number()` by 1 ulp                                                 | confirmed | [01](01-core-data-util.md#bugs)                                                       |
| 01/B2    | still open; partly fixed #2036 | Typed column `toArray({ array })` ignores the ctor and returns the source buffer             | [column-helpers.ts:35](../../packages/core/src/data/db/column-helpers.ts#L35)                                                                | BCIF conformation x/y/z alias the parsed Float64 columns: Float64 where Float32 was asked; a mutation would corrupt the source                    | confirmed | [01](01-core-data-util.md#bugs)                                                       |
| 01/B6    | still open   | `distinctColors` clustering reads `samples[]` where it means `sampleList[]`                  | [distinct.ts:144](../../packages/core/src/util/color/distinct.ts#L144)                                                                       | 'generate' palettes are biased (last `count` samples never assigned)                                                                              | confirmed | [01](01-core-data-util.md#bugs)                                                       |
| 03/B5    | still open   | Hash versions cover params only: tag and `dependsOn` edits never reach cells                 | [transform.ts:281](../../packages/core/src/state/transform.ts#L281)                                                                          | Hash-versioned updates (MVS) that change only tags or `dependsOn` have no effect on cells                                                         | confirmed | [03](03-core-state-task.md#bugs)                                                      |
| 03/B7    | still open   | A late abort discards a completed task's result; nobody disposes it (also 12/03/B7)          | [observable.ts:106](../../packages/core/src/task/execution/observable.ts#L106)                                                               | Abort after `task.f` returned: the built object (texture-mesh textures, a Download asset) leaks                                                   | confirmed | [03](03-core-state-task.md#bugs), [12](12-app-level.md#bugs)                          |
| 03/B2    | still open; partly fixed #2046 (no evaluation, replayed on a rebase; the flags reach the tree only) | Builder `updateState` never reaches `cell.state` and re-evaluates the subtree                | [builder.ts:357](../../packages/core/src/state/state/builder.ts#L357)                                                                        | Latent: no in-repo caller outside tests                                                                                                           | confirmed | [03](03-core-state-task.md#bugs)                                                      |

## v1 bugs now fixed

| Id     | PR    | Note                                                                                           |
| ------ | ----- | ---------------------------------------------------------------------------------------------- |
| 03/B1  | #1989 | flags carried through `flagLog` + `carryFlagWrites`; the fix introduced 03/B11                  |
| 04/B2  | #1989 | dynamic format properties for models 2..N                                                      |
| 06/B1  | #1989 | ASA `getValue` maps through `root.unitIndexMap`                                                |
| 05/B1  | #2022 | components share `a.root` (`remapModel(m, parent)`); the whole-unit sharing introduced 05/B9   |
| 07/N2  | #1995 | `setLayers(smooth = !merged)`, re-applied when merging changes                                 |
| 07/12/B1 | #1987 | isosurface ids kept on the unfilled space                                                    |
| 09/B1  | #2004 | intra `includeParent` filter reads `elements[members[i]]`                                      |
| 09/B4  | #2004 | bridge visual no longer rebuilds on every update                                               |
| 10/B1  | #2008 | the followed structure's sequence source gets its CIF frame                                    |
| 10/B2  | #2008 | `paintMarks` prunes removed representations                                                    |
| 10/B3  | #2008 | one interval per touched chain instance (10/B5 is a different gap there)                       |
| 11/B2  | #2008 | no snapshot within 30 s of one over the hard limit                                             |
| 14/N7  | #2014 | builder removal no longer writes into `state.tree` (`renewOwner`)                              |
| 02/B5  | #2025 | no coordinate read without a backbone atom (v2 bug; root of 02/P10)                     |
| 12/B2  | #2026 | v2 bug: the ESM entry's budget is on its start closure (5.5 MB; Prism now 5.16 MB)             |
| 11/B1  | #2029 | `For` has no `key` prop (a type error); the 49 sites pass `by=`                              |
| 11/B3  | #2029 | v2 bug: task pills follow their task by id (`lineTaskIds`, `lineTask`), made once             |
| 11/B4  | #2029 | v2 bug: History says the time of day, not "N min ago" written once                            |
| 11/B5  | #2029 | v2 bug: one Log bound (300) for every source, warnings and errors dropped last                |
| 08/B6  | #2028 | v2 bug (high): font atlas cache keyed by the font fields only, LRU of 4                        |
| 05/B9 | #2032 | `ofParent` also requires equal operators (structure.md §3) |
| 07/B2 | #2032 | new visuals get any non-identity repr transform (`isNew`) |
| 06/B12 | #2032 | the weak-H-bond refiner skips the intra check without intra contacts |
| 09/B7 | #2032 | relative interpolator; custom domain in sigma; direct volume and Prism `themeDomain` convert back |
| 09/B2 | #2032 | sequences by `byEntityKey`; rows without a sequence skipped |
| 04/B7 | #2032 | frames share the model's `IndexPairBonds` |
| 04/12/B3 | #2032 | periodic grids floor and wrap both neighbours |
| 06/B3 | #2032 | coarse elements get no radius or residue; applicable with an atomic unit only |
| 07/B3 | #2032 | coarse chains by asym id; residue/element levels skip coarse units |
| 01/B11 | #2032 | `o += 2 + xlen`; multi-member files still 01/B10 |
| 07/B1 | #2033 | a new volume or particle render object gets the representation state |
| 07/N1 | #2033 | segment mask reads zero outside the grid |
| 04/B3 | #2033 | ADPs mapped by atom id (`Map`); the provider is dynamic, per model |
| 09/B3 | #2033 | a set palette has at least 2 samples; past the last stop its color |
| 02/B1 | #2033 | `approxNearest` adds the nearest once (grid and structure lookups) |
| 06/B11 | #2033 | ligand negative charges from acid groups and charged atoms only |
| 04/B11 | #2033 | one bond per atom pair from either atom's CONECT; order the larger count |
| 04/B5 | #2033 | `'trip'` |
| 01/B12 | #2033 | a failed read of a held asset releases its hold |
| 03/B12 | #2033 | `updateTagged` sets params and tags both |
| 04/B6 | #2033 | `includes('Instances')` |
| 08/B7 | #2033 | the last dash gets a copy of the props |
| 03/B9 | #2034 | `subtree()` and `children()` skip refs without a cell |
| 03/B10 | #2034 | a stoppable walk returns the first match in pre-order |
| 05/B10 | #2034 | a copy with other elements gets its own indices; whole only for equal elements |
| 05/B12 | #2034 | units searched at `pivotR + maxRadius + radius` |
| 05/B6 | #2034 | `pick` pushes and pops the structure; `i % N === 0` |
| 06/B4 | #2034 | `indexX1 !== undefined` |
| 06/B7 | #2034 | legs and the donor–acceptor pair pass `validUnitPair` (`unitPairTest` for bridges undecided) |
| 07/B5 | #2034 | impostors rebuild on a size change with multiple or aromatic bonds or shortened cylinders |
| 08/B3 | #2034 | the weight is `v - floor(v)` (the cell is 1 or 0 wide) |
| 02/B2 | #2034 | chunks of ~2M visited voxels by the atom box |
| 10/B5 | #2034 | `wholeChains` uses the chain's atom span in its unit |
| 05/B2 | #2044 | coarse units read their own chain segments |
| 04/B4 | #2044 | `row + minSeqId`; the chain end is exclusive |
| 04/B12 | #2044 | bond atom ids mapped to rows; a bond to an unknown id dropped |
| 04/B13 | #2044 | a label with no set is skipped |
| 06/B6 | #2044 | the lookup radius from the largest radius assigned |
| 07/N3 | #2044 | the loci's instance, a set of labels per element |
| 08/B2 | #2044 | buckets and offsets over `maxId + 1` groups |
| 08/B4 | #2044 | vertex ids in a Uint32Array |
| 10/B6 | #2044 | a frame step starts the views without the clear and the replay |
| 12/B3 | #2044 | the last 1000 entries (`LogSize`) |
| 01/B16 | #2044 | a masked field's storage reads as its values (0 or '') |
| 03/B8 | #2045 | chunks sized from what was processed since the last update, a quarter per interval |
| 04/B14 | #2045 | masters kept per group name; a re-opened group gets its own |
| 09/B5 | #2045 | mesh and impostor `mustRecreate` ask the factory's support check |
| 09/B6 | #2045 | legs by `parentDisplay` on the end features; a shared leg drawn if any bridge draws it |
| 06/B9 | #2045 | an attach with the props of one in flight awaits it; computes itself if that one fails |
| 01/B3 | #2046 | `endI <= max` in the Interval branch |
| 01/B4 | #2046 | `set` of a present key replaces it and returns the old value |
| 01/B5 | #2046 | `is` tests `ctor`, `elementSize`, `chunks` (the padding half was #2007) |
| 01/B7 | #2046 | interval values tested; `projectValue` takes the segment holding the value |
| 01/B8 | #2046 | `start` clamped to `rowCount`, `end >= start` |
| 01/B9 | #2046 | typed arrays, Map and Set kept; other instances as before |
| 01/B13 | #2046 | doc: the gap-in-a tie order |
| 01/B15 | #2046 | `0x1f` mask, `fromCodePoint` for 4 bytes |
| 02/B3 | #2046 | extrema cleared when they cannot carry over |
| 02/B4 | #2046 | the state owns its centres |
| 02/12/B2 | #2046 | the stray line removed |
| 03/B3 | #2046 | `tags`/`dependsOn`/`state` actions; a replay adds to the newer tree |
| 03/B4 | #2046 | `setDependsOn` and `add` keep `tree.dependencies` |
| 05/B4 | #2047 | a `for` loop |
| 05/B5 | #2047 | the unit's own elements |
| 05/B8 | #2047 | unit operator then symmetry operator; the grid test in the close unit's frame |
| 05/B11 | #2047 | keyed by model id and operator |
| 06/B5 | #2047 | from key 0 |
| 06/B8 | #2047 | the slot-ordered `_energies` |
| 07/B6 | #2048 | from 0.475 to 0.5 of its own edge, as the mesh |
| 07/B7 | #2048 | `newGeometry || geometry`; shape positions before colors |
| 08/B5 | #2048 | `'3d'` |

Partly fixed: 03/B2 (#2046: no evaluation, replayed on a rebase; whether a builder's flags reach the cells is open). Dropped: 13/B1 (text mmCIF above ~4.8M atoms:
app paths decode to a `ChunkedBigString` since #2019; cost is now speed, 12/H5).
