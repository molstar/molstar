# Slice 01: core data and util

Scope: `packages/core/src/data` (int sets, db columns and tables, chunked arrays, buckets, iterators) and
`packages/core/src/util` (zip, memoize, LRU cache, object helpers, param definitions, colors). The code is mature and
rarely owns a hotspot outright. Its cost comes from how callers use it. Measured on `90d717f49` with Node 25.6 on an M4;
fixtures and stage numbers are in [00-profile.md](00-profile.md). Nothing here touches the frozen gl/Canvas3D code, but
`deepClone` (B9) is called from `gl/webgl/render-item.ts` and `gl/uniform.ts`.

## Hot paths

| Path | Where it shows up | Share |
| --- | --- | --- |
| `OrderedSet.getAt` / `indexOf` (Interval decode, binary search) per neighbour / per voxel | molecular-surface `obscured`, `projectTorus`; ASA; interactions; bonds-inter | big-bcif molsurf: getAt 881 ms + `Interval.getAt` 167 ms + `Tuple.is` 163 ms of 12.9 s; 7qpd molsurf `binarySearchRange` 6 %; ASA 6 %, interactions 13 %, bonds-inter 10 % |
| `ChunkedArray` `add3` / `allocateNext` / `_compact` | geometry builders | repr:interactions `_compact` 61 ms (about 85 % of it is `ret.set` memcpy of full chunks); ball-and-stick `allocateNext` 25 ms |
| `Table.view` → `arrayView`, `makeBuckets` via `sortAtomSite` | model stage | big-bcif `arrayView` 13.8 ms (15 %); big-cif `_makeBuckets` 11.3 ms self (12 %) |
| `ungzip` JS CRC-32 after native inflate | every `.gz` fetch / open | about 2 ms per MB uncompressed, on the main thread |

## Performance findings (confirmed and plausible)

| Id | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | `ungzip`: JS CRC after native inflate; native `gzip` stream 3.5x faster (and fixes B10, B11) | [zip.ts:231](../../packages/core/src/util/zip/zip.ts#L231) | time | medium | 190 of 275 ms per 100 MB; about 20-100 ms per real cif.gz/map.gz | S | confirmed | no |
| P2 | `OrderedSet.getAt`/`indexOf` re-dispatched in inner loops; no resolved accessor | [ordered-set.ts:37](../../packages/core/src/data/int/impl/ordered-set.ts#L37) | time | medium | hoist saves 1-1.5 ns of 5 ns per neighbour (about 25 % of `obscured`); 3-8 % of molsurf, ASA, interactions | M | confirmed | no |
| P3 | `subtract(Interval, SortedArray)`: binary search per element, overshoots to `b[end]` | [ordered-set.ts:299](../../packages/core/src/data/int/impl/ordered-set.ts#L299) | time | low | 4.8 ms vs 0.27 ms for 200k minus 100k (17x); a few ms per large deselect | S | confirmed | no |
| P4 | `makeBuckets`: `has` + `get` per element, second `getKey` pass | [buckets.ts:34](../../packages/core/src/data/util/buckets.ts#L34) | time | low | 4.3 ms → 1.4 ms per string pass at 186k atoms; about 6-10 ms of the 97 ms big-cif model stage | S | confirmed | no |
| P5 | `sortAtomSite`'s `Table.view` gathers every array column eagerly | [column.ts:485](../../packages/core/src/data/db/column.ts#L485) | time+memory | low | up to 13.8 ms and 10-20 MB transient at 186k sorted atoms, if only the used columns are materialised | M | plausible | no |
| P6 | `ChunkedArray._compact` copies a partial last chunk element by element | [chunked-array.ts:148](../../packages/core/src/data/util/chunked-array.ts#L148) | time | low | about 5 ms of 61 ms in big-bcif repr:interactions (the rest is the memcpy of full chunks) | S | confirmed | no |
| P7 | `SortedRanges` iterator: linear scan from range 0 per unit | [sorted-ranges.ts:92](../../packages/core/src/data/int/sorted-ranges.ts#L92) | time | low | 76 ms for 5000 units x 5000 ranges; 5-30 ms per polymer visual pass on models with more than 1000 chains | S | confirmed | no |
| P8 | `isParamEqual` on array values: `Object.keys` on both arrays | [param-definition.ts:755](../../packages/core/src/util/param-definition.ts#L755) | time+memory | low | 3.6 ms vs 0.23 ms at n=100k, 0.24 vs 0.02 ms at 10k; negligible at the sizes params hold today | S | confirmed | no |
| P9 | `distinctColors` uncached O(steps x count x samples) | [distinct.ts:117](../../packages/core/src/util/color/distinct.ts#L117) | time | low | 8 ms (75 colors) to 27-36 ms (250 colors) per theme build, only with the opt-in 'generate' palette | S | confirmed | no |

### P1: `ungzip` JS CRC

- Mechanism: the gzip header is parsed in JS and the raw blocks are inflated natively (`DecompressionStream('deflate-raw')`,
  [inflate.ts:169](../../packages/core/src/util/zip/inflate.ts#L169)). Then `crc(inflated)` runs byte by byte in JS on
  the main thread, and a mismatch only calls `console.error`.
- Evidence (`core-data-util/gz.ts`, 100 MB, level 1): `ungzip` takes 274-293 ms, of which the CRC is 192-202 ms.
  `new DecompressionStream('gzip')` takes 79-113 ms and validates the CRC itself. Callers:
  `plugin/core/src/state/transforms/data/fetch.ts:152` and `core/src/util/data-source.ts:146`.
- Fix: when `checkCompressionStreamSupport('gzip')` is true, pipe the whole file through `DecompressionStream('gzip')`.
  Write into the `isize` buffer while it fits, else concatenate the chunks. Keep the JS header path only as the
  fallback, fix B10 and B11 in it, and run the JS CRC there only.
- Risk: a corrupt file now throws instead of logging. Test truncated, multi-member and FEXTRA/FNAME files, and Node 22.

### P2: `OrderedSet` dispatch in inner loops

- Mechanism: `getAt` runs `I.is` (`typeof`) on every call; for an Interval it also decodes the IntTuple (Float64 store,
  Int32 load). `indexOf` on a SortedArray is a binary search. Callers run these per neighbour (`obscured`,
  `molecular-surface.ts:98`) and per voxel write (`projectTorus` `indexOf(indices, dp < 0 ? b : a)`, :311), where `a`
  and `b` are constant for the torus.
- Evidence (`verify/getat3.ts`, an `obscured`-shaped loop with 48 neighbours): getAt 4.9-5.5 ns vs hoisted 3.8-4.5 ns
  per neighbour, for both Interval and SortedArray sets. That matches the profile's roughly 9 % share in molsurf.
  The auditor's `getat2.ts` was too noisy to decide (equal after warm-up).
- Fix (this slice): an additive accessor, e.g. `OrderedSet.resolver(set) → { offset, array | undefined }`, and JSDoc
  saying `getAt`/`indexOf` are not for inner loops. Caller changes belong to their own slices: molsurf is in
  core-math (its `molsurf-opt.ts` already hoists `indexOf(a)`/`indexOf(b)` per torus); ASA, interactions and bonds
  are in model.
- Risk: none in core. Callers must keep resolved and unit-local indices apart.

### P3-P9 (low)

- **P3**: walk `v` from min to max with a pointer into `b` that starts at `findRange(b, min, max).start`, and stop at
  `max`. The overshoot (`last = b[min(end, li)]` with exclusive `end`) only triggers when `b` extends past `max(a)`.
  Unit-local loci and `polymer.ts:107` never do that, so it is latent: 68 ms for `[0..9] − [5, 2e7]`. The result is
  still correct because the typed array drops out-of-range writes.
- **P4**: add a prev-key fast path, use one `get` instead of `has` + `get`, and keep an Int32Array of bucket ordinals
  so the reorder pass needs no second `getKey`. Bench (`buckets.ts`): 4.3 → 1.4 ms. Only the entity and chain passes
  use string keys; `label_seq_id` is an int.
- **P5**: plausible, not measured. The saving depends on which permuted columns the model reads per atom. A lazy
  `arrayView` (gather on the first `toArray`/`__array`) has to keep `__array` truthy for `columnView` and
  `copyToArray`.
- **P6**: use `ret.set(lastChunk.subarray(0, currentIndex), offset)`. Line ticks in the big-bcif profile put about 8 % of
  `_compact` in the JS loop and about 85 % in `ret.set` of full chunks. The auditor's single over-allocated chunk case
  (3.1 → 0.3 ms) does not apply to the builders, which all pass `doNotResizeSingleton = true` and return `chunks[0]`.
  Geometric growth would cut `allocateNext` calls (25 ms in ball-and-stick); measure the peak memory first.
- **P7**: binary-search the first range with `max >= min(set)` and compare bounds directly, without an `Interval` per
  range. The existing top-of-call check already returns at once after a unit's last range, so the cost is the scan up
  to the unit's first range: O(U·R/2). Apply the same to `findSequenceKey` (`model/.../coarse-keys.ts:46`).
- **P8**: `Array.isArray` on both sides → `shallowEqualArrays`; `ArrayBuffer.isView` → compare length and elements.
- **P9**: an LRU memo keyed by `count` plus JSON of the props; swap-remove on an index array instead of `deepClone` +
  `splice`. Every built-in theme defaults to `'colors'`, and `particle-hierarchy` uses `clusteringStepCount: 0`, so
  the cost only appears when a user picks 'generate'.

## Bugs

| Id | Title | File | Severity | Verdict |
| --- | --- | --- | --- | --- |
| B11 | `ungzip` FEXTRA skip misses the 2-byte XLEN: every FEXTRA/BGZF gzip throws | [zip.ts:199](../../packages/core/src/util/zip/zip.ts#L199) | medium | confirmed (new) |
| B1 | `memoize1` never evicts; the `BasicData`-keyed caller leaks the source text per load | [memoize.ts:30](../../packages/core/src/util/memoize.ts#L30) | medium | confirmed |
| B10 | `ungzip` on multi-member gzip throws "inflates past its declared N bytes" | [zip.ts:226](../../packages/core/src/util/zip/zip.ts#L226) | low | confirmed (new) |
| B2 | Typed column `toArray({ array })` ignores the ctor and returns the source buffer | [column.ts:421](../../packages/core/src/data/db/column.ts#L421) | low | confirmed |
| B3 | `forEachSegment` emits the last segment twice for Interval sets | [ordered-set.ts:357](../../packages/core/src/data/int/impl/ordered-set.ts#L357) | low | confirmed |
| B6 | `distinctColors` clustering indexes `samples` instead of `sampleList` | [distinct.ts:150](../../packages/core/src/util/color/distinct.ts#L150) | low | confirmed |
| B9 | `deepClone` turns typed arrays, Map, Set and class instances into plain objects | [object.ts:90](../../packages/core/src/util/object.ts#L90) | low | confirmed (latent) |
| B4 | `LRUCache.set` appends a duplicate for an existing key; `get` is O(n) | [lru-cache.ts:45](../../packages/core/src/util/lru-cache.ts#L45) | low | confirmed (latent) |
| B5 | `ChunkedArray._compact` pads Array-backed output; `ChunkedArray.is` is always falsy | [chunked-array.ts:148](../../packages/core/src/data/util/chunked-array.ts#L148) | low | confirmed (latent) |
| B7 | `IntervalIterator` (start ≠ 0) and `Segmentation.projectValue` return wrong results | [interval-iterator.ts:32](../../packages/core/src/data/util/interval-iterator.ts#L32) | low | confirmed (latent) |
| B8 | `getArrayBounds` clamps `start` to `rowCount - 1`: an empty end window returns one row | [column-helpers.ts:9](../../packages/core/src/data/db/column-helpers.ts#L9) | low | confirmed (latent) |

- **B11**: `o += xlen` should be `o += 2 + xlen`. Failure: a gzip with FLG.FEXTRA (BGZF, some archivers) loses 2 bytes
  of the deflate stream and throws `TypeError` (`verify/gzextra.ts`); native `DecompressionStream('gzip')` decodes the
  same file. P1 fixes it.
- **B1**: `getUniqueComponentNames` (`model/.../basic/properties.ts:177`) is keyed by `BasicData`. Every load without
  `chem_comp` pins the CIF text: heap +4.1 MB per 7qpd load, retained through cache Map → `label_asym_id` → `str10` →
  the 1.8 MB string. `plugin/ui controls/parameters.tsx:962/964` adds an entry per edited color list. Fix: a `WeakMap`
  for object keys (branch on `typeof`), or compute in place.
- **B10**: `gzip -c a b > ab.gz` makes a valid file. `isize` comes from the last member, so the output overflows the
  buffer (`verify/gzmulti.ts`). P1 fixes it.
- **B2**: `getConformation` asks for Float32 x/y/z and gets the column's Float64Array: the BCIF frame buffer itself in
  identity order, or the `arrayView` copy when sorted (`coords.ts`: 3k7f and 1tup give Float64Array). Mutating
  `conformation.x` would corrupt the parsed column. Fix: copy into `params.array` when it differs from the source ctor.
  Model owners decide Float32 vs Float64 coordinates (B/atom vs exact-number tests).
- **B3**: `while (endI < _i …)` with an inclusive `_i`. `[0..3]` with one segment yields `[0, 0]`. The only caller
  (`plugin/ui sequence/polymer.ts:140`) deduplicates afterwards, so there is no visible effect today.
- **B6**: the closest-sample pass reads `samples[idx]` (an index into the spliced `sampleList`), and the assign pass
  covers `samples[0 .. n-count)`. The last `count` samples are never assigned, which biases the centroids. Fixing it
  changes 'generate' palettes slightly.
- **B9**: `'prototype' in source` is false for any instance: a Float32Array becomes `{0:…}`, a Map becomes `{}`, and a
  class instance loses its getters. A 1M-element Float32Array takes 54 ms. Today's only buffered uniform is a Vec2
  (Array), so this is latent. Fix in `util/object.ts` only: `ArrayBuffer.isView` → `slice()`; recurse only into plain
  objects.
- **B4**: the only caller (volume streaming) runs `queryData` inside state updates, which `State.updateQueue`
  serializes, and the camera and selection paths go through `SingleAsyncQueue`. The concurrent duplicate the auditor
  describes cannot happen today. A `Map<key, node>` makes `set` replace and `get` O(1).
- **B5, B7, B8**: no in-repo trigger (every ChunkedArray is typed; `IntervalIterator` and `projectValue` are only used
  in tests; no `toArray` caller passes `start`). Probes: Array-backed compact gives length 8 for 6 adds;
  `projectValue(5)` on `[0,10,20,30]` gives `(10,20)` and 25 gives `(30,0)`; `toArray({start:3,end:3})` gives `[3]`.

## Checked and dropped

- B4's "duplicate volume blocks under concurrent streaming queries": refuted, because every `queryData` call is
  serialized (see above). Only the latent API defect is kept.
- P6's "single over-allocated chunk" case and its "20-40 ms of the 61 ms" estimate: refuted, because builders pass
  `doNotResizeSingleton`, and the profile ticks are in the full-chunk memcpy.
- P7's "after its last range, scans on to the end": refuted, because the top-of-call `minAt(from) > max(set)` check
  returns at once. The quadratic part is the scan from range 0.
- P2 "1.5-2x per call" from `getat2.ts`: not reproducible after warm-up. The measured gain is about 25 % in a realistic
  loop body (`verify/getat3.ts`).
- P4's "initialise `prevKey` from `indices[start]`": harmless, because the first element always takes the new-bucket
  branch.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/core-data-util/`;
bundle each with `node ../build.mjs x.ts x.mjs && node x.mjs` (from `verify/` use `../../build.mjs`).

| Script | Settles |
| --- | --- |
| `gz.ts` | P1 ungzip vs JS CRC vs native gzip stream |
| `verify/gzmulti.ts`, `verify/gzextra.ts` | B10, B11 |
| `getat.ts`, `getat2.ts`, `verify/getat3.ts` | P2 (getat3 is the realistic `obscured`-shaped loop) |
| `sub.ts` | P3 subtract vs merge |
| `buckets.ts` | P4 |
| `chunked.ts` | P6 (plus the `_compact` line ticks from `../prof/big-bcif/repr_interactions.cpuprofile`) |
| `bugs.ts` | P3 overshoot, P7, B3, B5, B7 (`IntervalIterator`), B8 |
| `proj.ts` | B7 `projectValue` |
| `coords.ts` | B2 |
| `distinct.ts` | P9 |
| `verify/misc.ts` | B9 `deepClone`, P8 `isParamEqual`, B4 LRU duplicate |
