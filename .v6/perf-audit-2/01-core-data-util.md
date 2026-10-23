# Slice 01: core data, util and algorithms (audit v2)

Scope: `packages/core/src/data` (int sets, db columns and tables, chunked arrays, buckets, iterators),
`packages/core/src/util` (zip, number parser, utf8, memoize, LRU cache, object helpers, param definitions, colors,
assets, html entities) and, new since v1, `packages/core/src/algorithms` (`sequence/needleman-wunsch.ts`). Run at
`0944823d4` with Node 25.6 on an M4. Fixtures and stage numbers are in [00-profile.md](00-profile.md). The slice still
rarely owns a hotspot outright: its cost comes from the model stage and the geometry builders. Nothing here touches the
frozen gl/Canvas3D code.

## Hot paths now

| Path                                                                 | Where it shows up                              | Share                                                                                                                                                    |
| -------------------------------------------------------------------- | ---------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Table.view` → `arrayView` (`sortAtomSite`, sorted models)           | model stage                                    | flat-192 (1.12M atoms): 97 ms self, 17 % of the stage (about 75 ms in the gather loop, about 20 ms on `new ctor(n)`); 13 had 16 % of flat-864's 5 s |
| `makeBuckets` (entity, chain, seq passes)                            | model stage                                    | flat-192: `_makeBuckets` 39 ms + `makeBuckets` 15 ms (inlined) self; 7qpd 0.5 ms (5 %)                                                                   |
| `ChunkedArray.add3` / `add` / `_compact` / `create`                  | geometry builders                              | flat-192 ball-and-stick: `add3` 190 + `add` 150 ms self (8 % of 4.4 s); interactions: `_compact` 108, `add3` 76, `create` 67, `add` 57 ms                |
| `SortedRanges` iterator (`firstIntersectionIndexFrom` from range 0)  | structure stage (`getAtomicPolymerElements`)   | flat-192: about 6 ms of a 40 ms stage (`firstIntersectionIndexFrom` 1.7, `Interval.ofRange` 1.7, `areIntersecting*` 2.6 ms)                               |
| `OrderedSet.getAt` via core-math `getBoundary`                       | lookup+boundary                                | 7qpd: `getAt` 1.7 of about 7.6 ms (22-25 %); the molsurf and GridLookup3D callers were converted in #1988 and #1992                                       |
| `ungzip` JS CRC-32 after native inflate                              | every `.gz` fetch or open                      | about 2 ms per MB uncompressed, main thread                                                                                                              |

Changed since v1 in these paths: #2007 (`ChunkedArray._compact` copies only the used part of the last chunk; new
`compactTrimmed` used by the mesh, lines, spheres and cylinders builders), #2017 (shared Gotoh DP in
`algorithms/sequence/needleman-wunsch.ts`; `sortAtomSite` reads in-order models in place), #2019 (`AssetManager`
interface, `DefaultAssetManager`), #2021 (`Clip` tag), #2008 (`html-entities`, Markdown decode, `escapeHtml`), `unzip`'s
`only`. `ungzip`, `OrderedSet`, `SortedRanges`, `buckets`, `column*`, `memoize`, `lru-cache`, `object`, `distinct` and
`param-definition` are unchanged.

## v1 → now

Fixed 1, partly fixed 2, still open 17, regressed 0, dropped 0 (20 v1 items: P1-P9, B1-B11).

| v1 id | Title                                                                        | Status       | PR    | Note                                                                                                                                         |
| ----- | ---------------------------------------------------------------------------- | ------------ | ----- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| P1    | `ungzip`: JS CRC after native inflate                                        | still open   |       | zip.ts:199-250 unchanged; `crc()` over the whole output at :245, a mismatch only logs                                                       |
| P2    | `OrderedSet.getAt`/`indexOf` re-dispatched in inner loops                    | still open   |       | no resolver API; the big callers (molsurf `obscured`, GridLookup3D) hoisted it themselves in #1988/#1992; `getBoundary` still 22-25 % of 7qpd |
| P3    | `subtract(Interval, SortedArray)`: binary search per element, overshoot      | still open   |       | `subtractIS` ordered-set.ts:299-329 unchanged                                                                                               |
| P4    | `makeBuckets`: `has` + `get` per element, second `getKey` pass               | still open   |       | buckets.ts:32-43, :75-79 unchanged; 39 + 15 ms self on flat-192                                                                             |
| P5    | `sortAtomSite`'s `Table.view` gathers every array column eagerly             | partly fixed | #2017 | in-order models read typed columns in place (sort.ts:47-59); sorted models still `Table.view` (sort.ts:63)                                   |
| P6    | `ChunkedArray._compact` copies a partial last chunk element by element       | fixed        | #2007 | `ret.set(lastChunk.subarray(0, currentIndex))` at chunked-array.ts:155; the Array path loops only to `currentIndex`                          |
| P7    | `SortedRanges` iterator: linear scan from range 0 per unit                   | still open   |       | sorted-ranges.ts:92-105 unchanged                                                                                                           |
| P8    | `isParamEqual` on array values via `Object.keys`                             | still open   |       | param-definition.ts:754-755 unchanged                                                                                                       |
| P9    | `distinctColors` uncached                                                    | still open   |       | distinct.ts:117-170 unchanged                                                                                                               |
| B1    | `memoize1` never evicts; `BasicData`-keyed caller leaks the source per load  | still open   |       | memoize.ts:29-37 a Map; properties.ts:189 keyed by `BasicData`; plugin/ui parameters.tsx:962/964 (keyed by color arrays) are still there   |
| B2    | Typed column `toArray({ array })` ignores the ctor, aliases the source       | still open   |       | column.ts:421 → `typedArrayWindow` column-helpers.ts:35-40 unchanged                                                                         |
| B3    | `forEachSegment` emits the last segment twice for Interval sets              | still open   |       | ordered-set.ts:357 `endI < _i`, `_i` inclusive                                                                                              |
| B4    | `LRUCache.set` appends a duplicate for an existing key                       | still open   |       | lru-cache.ts:45-55 unchanged; only caller volume-streaming behavior.ts:294                                                                   |
| B5    | `_compact` pads Array-backed output; `ChunkedArray.is` always falsy          | partly fixed | #2007 | padding fixed (chunked-array.ts:160); `is` (:34-36) still tests `creator`/`chunkSize`; no callers                                            |
| B6    | `distinctColors` clustering indexes `samples` instead of `sampleList`        | still open   |       | distinct.ts:144, :148-151                                                                                                                   |
| B7    | `IntervalIterator` (start ≠ 0) and `Segmentation.projectValue` wrong         | still open   |       | interval-iterator.ts:21-35, segmentation.ts:52 unchanged; tests only                                                                        |
| B8    | `getArrayBounds` clamps `start` to `rowCount - 1`                            | still open   |       | column-helpers.ts:10                                                                                                                        |
| B9    | `deepClone` turns typed arrays, Map, Set, instances into plain objects       | still open   |       | object.ts:90                                                                                                                                |
| B10   | `ungzip` on multi-member gzip throws                                         | still open   |       | zip.ts:240-241 `new Uint8Array(isize)` of the last member                                                                                   |
| B11   | `ungzip` FEXTRA skip misses the 2-byte XLEN                                  | still open   |       | zip.ts:213 `o += xlen`                                                                                                                      |

## Performance findings (confirmed and plausible, open)

| Rank | Id  | vs v1        | Title                                                                                     | File                                                                                    | Kind        | Impact | Estimate                                                                       | Effort | Verdict   | Frozen |
| ---- | --- | ------------ | ----------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- | ----------- | ------ | ------------------------------------------------------------------------------ | ------ | --------- | ------ |
| 1    | P10 | new; **fixed** #2032 | `arrayView` gather: one polymorphic keyed load/store for every column kind                | [column.ts:485](../../packages/core/src/data/db/column.ts#L485)                         | time        | medium | -40 to -50 ms of 97 ms on flat-192 model; about -0.2 to -0.4 s on flat-864      | S      | confirmed | no     |
| 2    | P1  | still open   | `ungzip` runs a JS CRC-32 over the output after native inflate                            | [zip.ts:245](../../packages/core/src/util/zip/zip.ts#L245)                              | time        | medium | 190 of 275 ms per 100 MB; 20-100 ms per real cif.gz/map.gz                      | S      | confirmed | no     |
| 3    | P12 | new          | `ChunkedArray.add`/`add3`/`add4` stores are 4-way polymorphic across builders             | [chunked-array.ts:59](../../packages/core/src/data/util/chunked-array.ts#L59)           | time        | low    | up to -100 ms of 340 ms on flat-192 ball-and-stick; -40 ms interactions         | M      | plausible | no     |
| 4    | P4  | still open   | `makeBuckets`: `has` + `get` per element, second `getKey` pass                            | [buckets.ts:34](../../packages/core/src/data/util/buckets.ts#L34)                       | time        | low    | -15 to -25 ms on flat-192 model; about -80 ms at flat-864                       | S      | confirmed | no     |
| 5    | P11 | new          | Needleman-Wunsch: O(n·m) traceback unbounded in the model path; matrix kind polymorphic   | [needleman-wunsch.ts:80](../../packages/core/src/algorithms/sequence/needleman-wunsch.ts#L80) | time+memory | low    | +27 % from mixed matrix kinds; a 5k chain is 25 MB and 135 ms, 20k is 400 MB and 2.1 s | M      | confirmed | no     |
| 6    | P5  | partly fixed | Sorted models: `Table.view` copies every array column of `atom_site` eagerly               | [sort.ts:63](../../packages/model/src/formats/structure/basic/sort.ts#L63)              | time+memory | low    | beyond P10: the transient copies of rarely read columns (about 9 MB each at 1.1M atoms) | M      | plausible | no     |
| 7    | P2  | still open   | `OrderedSet.getAt`/`indexOf` re-dispatch and decode per call                              | [ordered-set.ts:37](../../packages/core/src/data/int/impl/ordered-set.ts#L37)           | time        | low    | up to 25 % of small-entry boundary time (1.7 ms on 7qpd)                        | M      | confirmed | no     |
| 8    | P7  | still open   | `SortedRanges.firstIntersectionIndexFrom` scans linearly from range 0 per unit            | [sorted-ranges.ts:92](../../packages/core/src/data/int/sorted-ranges.ts#L92)            | time        | low    | about 6 ms per polymer pass at 2.5k chains (flat-192); O(U·R/2)                 | S      | confirmed | no     |
| 9    | P3  | still open   | `subtractIS` does a binary search per element of the interval                            | [ordered-set.ts:326](../../packages/core/src/data/int/impl/ordered-set.ts#L326)         | time        | low    | 4.8 vs 0.27 ms for 200k minus 100k; a few ms per large deselect                 | S      | confirmed | no     |
| 10   | P9  | still open   | `distinctColors` uncached O(steps × count × samples)                                      | [distinct.ts:117](../../packages/core/src/util/color/distinct.ts#L117)                  | time        | low    | 8-36 ms per theme build, opt-in 'generate' only                                 | S      | confirmed | no     |
| 11   | P8  | still open   | `isParamEqual` compares array and typed-array values via `Object.keys`                    | [param-definition.ts:755](../../packages/core/src/util/param-definition.ts#L755)        | time+memory | low    | negligible at today's param sizes                                               | S      | confirmed | no     |

### P10: `arrayView` is polymorphic over column kinds

- Mechanism: `arrayView` does `new (array.constructor)(n)` and `ret[i] = array[map[i]]` in one function for every
  column `Table.view` permutes: Float32Array, Float64Array, Int32Array and, for BCIF string fields (`__array` is the
  decoded `string[]`, so `columnView`'s test at :481 passes), plain Arrays. The keyed load and store see every kind.
- Evidence: `aview.ts` (22 columns × 1.1M rows): `Column.view` 86 ms in round 1, 441-458 ms in rounds 2-4; per-kind
  monomorphic copies 44-46 ms. `v/aview3.ts` splits it: typed columns alone stay at 1.2-2.8 ms/col; once string Arrays
  share the site, typed columns cost 25-34 ms/col (15-20×), strings 6 ms. The profiled flat-192 run sat between the two
  (97 ms for about 20 permuted columns, about 4 ms/col), so the microbench cliff is the worst case, not what the pipeline
  measured. Line ticks: about 75 ms in the loop, about 20 ms on `new ctor(n)`.
- Fix: dispatch once per column to per-kind gathers (Float32, Float64, Int32, Uint8/Int8, `Array.isArray` → `new
  Array(n)`, generic), each a separate source function so each keeps its own feedback; a `gatherInto` helper in
  `data/util/array.ts`.
- Risk: low; keep `__array` truthy and `valueKind` as now (`columnView`, `copyToArray`, the hierarchy fast paths).

### P1: `ungzip` JS CRC

- Mechanism: the gzip header is parsed in JS, blocks inflated by `DecompressionStream('deflate-raw')`, then `crc()`
  walks every byte on the main thread; a mismatch only logs.
- Evidence: v1 `gz.ts`, 100 MB: 274-293 ms, CRC 192-202 ms; native `DecompressionStream('gzip')` 79-113 ms with the CRC
  checked. Code unchanged; callers `plugin/core/.../data/fetch.ts:152`, `core/util/data-source.ts:146`.
- Fix: with `checkCompressionStreamSupport('gzip')`, pipe the whole file through `DecompressionStream('gzip')`, writing
  into the `isize` buffer while it fits, else concatenating. Keep the JS path only as the fallback, with B10 and B11
  fixed. Fixes B10 and B11 on the native path.
- Risk: a corrupt file throws instead of logging; test truncated, multi-member, FEXTRA and FNAME files on Node 22.

### P12: `ChunkedArray` stores polymorphic (low)

- Mechanism: one shared `add3`/`add`/`add4` writes into `currentChunk` for every builder. Non-test ctors: Float32Array 28,
  Uint32Array 10, Int32Array 3 (OBJ parser only), Uint8Array 2 (interaction features); Array only in tests. A session
  sees 2-4 kinds, so the keyed stores are polymorphic, not megamorphic.
- Evidence: `chunk.ts`, 3M `add3` into Float32: 9.0-9.2 ms monomorphic, 13.1-13.5 ms after Uint32, Int32 and Uint8
  (+45 %); with no OBJ loaded, three kinds, likely a little less. The share in real builders is extrapolated from the profile (`add3` 190 + `add` 150 ms in flat-192
  ball-and-stick), not A/B measured, hence plausible.
- Fix: per-kind entry points as separate source functions (`add3F32`, `addU32`, ...). Closures made by one factory
  share their feedback cell in V8, so a factory returning a closure per ctor does not split the feedback. Separately,
  look at `create`'s 67 ms in interactions (zero-filled estimate-sized chunks, trimmed later by `compactTrimmed`).
- Risk: API churn in the graphics builders; no behaviour change.

### P4, P11, P5, P2, P7, P3, P9, P8 (low)

- **P4**: line ticks put about 57 of 255 `_makeBuckets` ticks on `has` + `get` and about 50 on the reorder pass. Fix: a
  `prevKey === key` fast path on the cached bucket, one `get` instead of `has` + `get`, and an Int32Array of bucket
  ordinals from pass 1 so pass 2 needs no `getKey` or Map lookup. Same buckets.
- **P11**: `globalAlign` allocates a (n+1)(m+1) Uint8Array traceback. Seq* caps at 2^26 cells (`checkCells`); the
  model's `align` (alignment.ts:73) does not. Its callers are PDB SEQRES mapping (atom-site.ts:110, per chain) and
  sequence superposition (sequence.ts:60). `nw.ts`: 1k 9 ms, 3k 45 ms, 10k 527 ms, 20k 2,107 ms (5.0-5.3 ns/cell,
  567 MB RSS at 20k); after the Int8 matrix (Seq*) and the no-matrix path ran, 6.6-6.9 ns/cell (+27 %). `nw2.ts`'s
  query profile gave no inner-loop gain. `nwfuzz.ts`: 24,000 brute-force cases, 0 mismatches. Real PDB chains stay
  under about 5k residues (about 135 ms, 25 MB), so this is a guard, not a hotspot. Fix: one matrix kind for both
  callers (Float64 or Int32), a cell cap in the model path, and a band or Hirschberg fallback above a budget.
- **P5**: plausible, low. Most permuted columns are read by the hierarchy and kept as the model's atom columns
  (atomic.ts:67-152); `Cartn_*.toArray({ array: Float32Array })` returns the `arrayView` copy itself (B2). So a lazy
  view saves only the transient copies of rarely read columns (`type_symbol` is mapped to an ElementSymbol array right
  away), and lazy `value` reads add a `map[row]` indirection to every later read. The time half is P10.
- **P2**: add an additive `OrderedSet.resolver(set)` (or `forEachIndexed`) and JSDoc saying `getAt`/`indexOf` are not
  for inner loops; convert `getBoundary` (core-math boundary.ts:44-60) and the remaining per-element callers.
- **P7**: binary-search the first range with `max >= min(set)`, compare bounds against `findPredecessorIndex` without an
  `Interval` per range; same in model `coarse-keys.ts` `findSequenceKey`.
- **P3**: two-pointer walk from `findRange(b, min, max).start`, stopping at `max` (the `last = b[min(end, li)]`
  overshoot is latent; typed-array writes past the end are dropped, so the result is right).
- **P9**: LRU memo keyed by count + JSON of props; swap-remove on an index array (fix B6 with it).
- **P8**: `Array.isArray` on both → `shallowEqualArrays`; `ArrayBuffer.isView` → length and elements.

## Bugs

| Rank | Id  | vs v1        | Severity | Title                                                                                       | File                                                                                    | Verdict   |
| ---- | --- | ------------ | -------- | ------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- | --------- |
| 12   | B1  | still open   | medium   | `memoize1` never evicts; `getUniqueComponentNames` keyed by `BasicData` pins each source    | [properties.ts:189](../../packages/model/src/formats/structure/basic/properties.ts#L189), [memoize.ts:29](../../packages/core/src/util/memoize.ts#L29) | confirmed |
| 13   | B11 | **fixed** #2032    | medium   | `ungzip` FEXTRA skip omits the 2-byte XLEN: every FEXTRA/BGZF `.gz` fails                   | [zip.ts:213](../../packages/core/src/util/zip/zip.ts#L213)                              | confirmed |
| 14   | B12 | new; **fixed** #2033 | low      | `AssetManager.resolve` keeps a reference when reading the asset fails                       | [assets.ts:225](../../packages/core/src/util/assets.ts#L225)                            | confirmed |
| 15   | B10 | still open   | low      | `ungzip` trusts the last member's ISIZE: multi-member and >4 GiB files throw                | [zip.ts:241](../../packages/core/src/util/zip/zip.ts#L241)                              | confirmed |
| 16   | B14 | new          | low      | `number-parser`: `parseInt` wraps above 2^31-1; `parseFloat` is not correctly rounded       | [number-parser.ts:39](../../packages/core/src/util/number-parser.ts#L39)                | confirmed |
| 17   | B16 | new          | low      | `Column.view`/`Column.window` of a masked typed column return the raw slot, not the default | [column.ts:481](../../packages/core/src/data/db/column.ts#L481)                         | confirmed |
| 18   | B2  | still open   | low      | Typed column `toArray({ array })` ignores the ctor and returns the source buffer            | [column-helpers.ts:35](../../packages/core/src/data/db/column-helpers.ts#L35)           | confirmed |
| 19   | B6  | still open   | low      | `distinctColors` clustering reads `samples[]` where it means `sampleList[]`                 | [distinct.ts:144](../../packages/core/src/util/color/distinct.ts#L144)                  | confirmed |
| 20   | B3  | still open   | low      | `OrderedSet.forEachSegment` emits the last segment twice for Interval sets                  | [ordered-set.ts:357](../../packages/core/src/data/int/impl/ordered-set.ts#L357)         | confirmed |
| 21   | B13 | new          | low      | `globalAlign` documents a gap-tie order the D state does not follow                        | [needleman-wunsch.ts:137](../../packages/core/src/algorithms/sequence/needleman-wunsch.ts#L137) | confirmed |
| 22   | B15 | new          | low      | `_utf8Read` fallback decodes 2-byte and 4-byte UTF-8 wrongly (latent)                       | [utf8.ts:74](../../packages/core/src/util/utf8.ts#L74)                                  | confirmed |
| 23   | B9  | still open   | low      | `deepClone` turns typed arrays, Map, Set and instances into plain objects (latent)          | [object.ts:90](../../packages/core/src/util/object.ts#L90)                              | confirmed |
| 24   | B4  | still open   | low      | `LRUCache.set` appends a duplicate for an existing key (latent)                             | [lru-cache.ts:45](../../packages/core/src/util/lru-cache.ts#L45)                        | confirmed |
| 25   | B5  | partly fixed | low      | `ChunkedArray.is` is always falsy (latent)                                                  | [chunked-array.ts:34](../../packages/core/src/data/util/chunked-array.ts#L34)           | confirmed |
| 26   | B7  | still open   | low      | `IntervalIterator` (start ≠ 0) and `Segmentation.projectValue` wrong (latent)               | [interval-iterator.ts:21](../../packages/core/src/data/util/interval-iterator.ts#L21)   | confirmed |
| 27   | B8  | still open   | low      | `getArrayBounds` clamps `start` to `rowCount - 1` (latent)                                  | [column-helpers.ts:10](../../packages/core/src/data/db/column-helpers.ts#L10)           | confirmed |

- **B1**: `memoize1` caches in a Map forever. `getUniqueComponentNames` (used at properties.ts:136, :181) is keyed by
  `BasicData`. Failure: every mmCIF load without `chem_comp` pins atom_site → column closures → the source text or
  frame (v1: +2.7-4.1 MB per 7qpd load, retained). plugin/ui `parameters.tsx:962/964` adds an entry per color array
  it renders (v1 said these were gone: they are not). Fix: a WeakMap for object keys in `memoize1`, or compute in place.
- **B11**: `o += xlen` should be `o += 2 + xlen`. Failure: any FLG.FEXTRA gzip (BGZF, some archivers) starts the
  deflate stream 2 bytes early and throws (v1 `verify/gzextra.ts`); native `DecompressionStream('gzip')` decodes it.
- **B12**: the already-held paths (:225 URL, :263 file), the uncached POST path (:237) and the stored-File path (:270)
  take a reference before `await readFromFile(...)`, with no catch; only the byte-cache path (:252-257) releases on
  failure. Failure (`v/assets-leak.ts`): a File asset whose body is not JSON, resolved as `'json'`, throws and stays
  held (`assets.length` 1). Held once as `'string'`, then failed as `'json'`, then the good wrapper disposed: the entry
  stays (refCount stuck at 1). Real triggers: a GraphQL POST answered by an HTML error page (graphql-client.ts:69), a
  dropped FEXTRA `.gz` (B11) or a corrupt file, an aborted read. Fix: try/catch → `this.release(asset)` and rethrow
  around each such read.
- **B10**: `buf = new Uint8Array(isize)` with `isize` from the last 4 bytes (last member only, modulo 2^32).
  Failure: `gzip -c a b > ab.gz` throws "inflates past its declared N bytes" or trailing-junk (v1 `verify/gzmulti.ts`);
  a map.gz over 4 GiB uncompressed gets `isize = size mod 2^32` and throws the same way. Fix: P1's native path with
  chunk concatenation; in the JS fallback, loop over members.
- **B14**: `parseInt` accumulates `(10 * ret + c) | 0`. Failure (`pf.ts`): `'3000000000'` → -1294967296,
  `'2147483648'` → -2147483648, `'-2147483649'` → 2147483647, silently. `parseFloat` computes `ret + point / div` and
  `main * 10^exp` (two roundings): 149 of 200,000 random 3-decimal values differ from `Number()` by 1 ulp
  (`'1.805'` → 1.8050000000000002) and 50,864 of 200,000 scientific values (`'8.1110e14'` → 811100000000000.1).
  Int fields above 2^31 are rare in mmCIF, so this is a precision and robustness defect. Fix: accumulate in a double
  and fall back to `Number()` past 9 digits; for floats, Clinger's fast path (integer mantissa < 2^53, |k| ≤ 22, one
  exact power-of-ten multiply or divide), else `Number(substring)`. Risk: last-bit changes in exact-number snapshots.
- **B16**: `columnView` (:481) and `windowTyped` (:442) read `__array` directly, but a masked BCIF numeric field's
  `value` returns 0 for masked rows, not the stored slot (binary/field.ts:32-46). Failure (`v/mask.ts`): data
  `[1, 7, 3, 4]`, row 1 masked: the source `value(1)` is 0, `Column.view(c, [3,1,2,0]).value(1)` and
  `Column.window(c, 1, 3).value(0)` are 7. Mol*'s own encoder leaves masked numeric slots 0, so only BCIF from an
  encoder that fills them otherwise shows it (`label_seq_id` of HETATM, `pdbx_formal_charge`) in sorted or windowed
  `atom_site`s. Fix: when the source column has a mask, gather with `typedArrayWindowMasked`-style defaults, or let a
  field expose whether `value` is a plain array read.
- **B2**: `typedArrayWindow` ignores `params.array` and returns `data` for the full range. Failure: `getConformation`
  asks for Float32 x/y/z and gets the Float64 column buffer (the BCIF frame or the `arrayView` copy); mutating it would
  corrupt the parsed column. Fix: copy into `params.array` when it differs from `data.constructor`.
- **B6**: zones get `samples[idx]` with `idx` into the spliced `sampleList`, and the assign loop covers
  `samples[0 .. sampleList.length)`: the last `count` samples are never assigned, so 'generate' palettes are biased.
- **B3**: `[0..3]` with one segment yields `[0, 0]`; the only caller (plugin/ui `sequence/polymer.ts`) deduplicates,
  so it is masked today. Fix: `endI <= _i`.
- **B13**: the doc (:39-41) says `'open'` prefers a pair, then extending, then the other gap. The I state does (:129).
  The D state under `'open'` takes the extension only when strictly better than both opens (:137), so its order is
  pair, other gap, extend. This matches the old Seq* argmax the results are pinned to, so only the doc is wrong.
  Failure: a new caller that relies on the documented ties gets a different one of the equal-score alignments. Fix the
  doc.
- **B15**: the two-byte branch masks the lead with 0x0f instead of 0x1f (U+0400-U+07FF decode wrong: Cyrillic, Hebrew,
  Arabic), and the four-byte branch passes a code point above U+FFFF to `String.fromCharCode`. Used only without
  `TextDecoder` (:104-108), so latent. Fix: 0x1f with `String.fromCharCode` (the 1024-entry table is too small for
  0x7ff), `String.fromCodePoint` for four bytes, or drop the fallback.
- **B9, B4, B5, B7, B8**: latent, unchanged since v1. `deepClone` copies keys of any instance (a Float32Array becomes
  `{0: ...}`); `LRUCache.set('a')` twice keeps two entries (streaming is serialized); `ChunkedArray.is` tests
  `creator`/`chunkSize`, which no ChunkedArray has; `projectValue(5)` on `[0,10,20,30]` gives (10, 20);
  `toArray({ start: 3, end: 3 })` gives one row.

## Checked and dropped

- P10's "4-9× cliff on later loads": the cliff is real in the microbench (typed columns 15-20× slower once string
  Arrays share the site) but the profiled flat-192 run did not hit it (about 4 ms/col). Kept as the worst case only.
- P5's "tens of MB retained if rarely read columns stay lazy": most permuted columns are the model's own atom columns,
  and the coordinates' `toArray` returns the copy. The memory saving is transient copies of a few columns.
- P12's "closure-specialized adder per ctor" fix: closures from one function literal share feedback in V8; the fix
  needs separate source functions.
- P11's titin-size estimate (35k × 35k, 1.2 GB, 6.6 s): extrapolated, and no PDB chain is that long; the measured
  sizes stand.
- v1 B1 evidence "the plugin/ui parameters.tsx callers are gone": wrong, both are at parameters.tsx:962/964.
- P1 "low" (auditor): kept at medium as in v1; the code and the 2 ms/MB main-thread cost are unchanged.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/core-data-util/`;
build each with `node build.mjs $PWD/x.ts $PWD/x.mjs && node x.mjs`.

| Script                     | Settles                                                                 |
| -------------------------- | ----------------------------------------------------------------------- |
| `aview.ts`, `aview2.ts`    | P10 `Column.view` vs per-kind gathers, mixed and typed-only             |
| `v/aview3.ts`              | P10 per-kind breakdown (`all`, `typed`, `stronly`)                      |
| `chunk.ts`                 | P12 `add3` monomorphic vs after four ctors (`mono`, `poly`)             |
| `nw.ts`, `nw2.ts`          | P11 cells/time/RSS by size (`mono`, `poly`); query-profile variant      |
| `nwfuzz.ts`                | `globalAlign` vs brute force, all tie and gap modes (24,000 cases)      |
| `pf.ts`                    | B14 `parseInt` wrap, `parseFloat` vs `Number()`                         |
| `v/assets-leak.ts`         | B12 failed reads stay held                                              |
| `v/mask.ts`                | B16 masked column through `Column.view` and `Column.window`             |
| `../profile/prof-flat-192` | line ticks for `arrayView`, `_makeBuckets`; structure-stage SortedRanges |

v1 scripts (P1 `gz.ts`, B10/B11 `verify/gzmulti.ts`, `verify/gzextra.ts`, P3 `sub.ts`, P2 `verify/getat3.ts`, B2
`coords.ts`, P9 `distinct.ts`, B9/P8/B4 `verify/misc.ts`) are listed in [../perf-audit/01-core-data-util.md](../perf-audit/01-core-data-util.md).
