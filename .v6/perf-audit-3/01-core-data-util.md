# Slice 01: core data, util and algorithms (v3)

Scope: `packages/core/src/data` (int sets, db columns and tables, chunked arrays, buckets, iterators),
`packages/core/src/util` (zip, number parser, utf8, sniffable, data-source, assets, memoize, LRU cache, object helpers,
param definitions, colors and the new value palettes, text measure) and `packages/core/src/algorithms`
(`sequence/needleman-wunsch.ts`). HEAD `16277d822`, Node 26.12 on the M4. **Machine load 26–33** during the runs (other
auditors): wall times are 2–3× v2's and carry ±50 %; quantitative claims below use min-of-N, same-process A/B,
allocation (`total_allocated_bytes`) or `arrayBuffers`. A/B base: `git archive 0944823d4 packages/core` bundled by an
esbuild plugin that maps `@molstar/core/*` to it (`build-base.mjs`). Nothing here touches frozen gl/Canvas3D code.

PRs since v2 in these paths: #2032 (`arrayPickIndices` per kind, FEXTRA), #2033 (asset holds on failed reads),
#2035 (`fetchData` by `fetch`), #2036 (`Column.storage`, dictionary columns, `parse*Bytes`, utf8 short-ASCII path),
#2039 (assets hold bytes, `decompress` by magic, `Sniffable`, `asciiIndexOf`, `isGzip`), #2040 (`ReadonlyUint8Array`,
utf8/number-parser byte variants), #2044 (`ungzip` by native `DecompressionStream('gzip')`), #2045 (`makeBuckets`),
#2046 (14 latent bugs, 10 perf: ordered-set, sorted-ranges, segmentation, interval-iterator, LRU, distinct cache,
deepClone, isParamEqual, utf8, async-queue), #2047 (`arrayStats`), #2043 (new `value-palette`, `split-color`),
#2030 (new `text-measure`).

## Hot paths now

| Path                                                                                          | Fixture, stage               | Cost (loaded machine, profile)                                                                                        |
| --------------------------------------------------------------------------------------------- | ---------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `Column.view` → `columnView` (incl. lazy BCIF decode on `storage()`), `arrayPickIndices` per kind | flat-192 model (1.33 s)       | `columnView` 214 ms total (16 %), of which `arrayPickIndices` 75 ms (`pickUint8` 39, `pickInt32` 27, `pickUint16` 9 self) |
| `makeBuckets` / `_makeBuckets`                                                                | flat-192 model               | 126 ms total (10 %), `_makeBuckets` 67 ms self                                                                        |
| `typedArrayWindowMasked` (masked BCIF storage copies), `Segmentation.create`                  | flat-192 model               | 28 + 28 ms self                                                                                                       |
| `asciiIndexOf` over the bytes (CIF variant guess, 5 needles absent in a plain mmCIF)           | 124 MB flat-192.cif, `auto`  | 603 ms min-of-7 (vs 268 ms `string.indexOf` on the same text)                                                         |
| `SortedRanges` iterator + per-residue `areIntersecting` in `getAtomicPolymerElements`         | flat-192 structure (179 ms)  | 94 ms total: `binarySearchPredIndexRange` 38, `areIntersecting*` 37, `Interval.ofRange` 21 ms                          |
| `OrderedSet.getAt` in core-math `getBoundary`                                                 | flat-192 lookup+boundary     | `getAt2` 30 ms total of 523 ms `getBoundary`                                                                          |
| `ungzip` (native gzip now)                                                                    | 50 MB synthetic `.gz`        | 200 ms HEAD vs 1,520 ms v2 base, same process load                                                                    |

## v2 → now

Fixed 16, partly fixed 1, still open 10, regressed 0, dropped 0 (27 v2 items: P1–P5, P7–P12, B1–B16).

| v2 id | Title                                                              | Status       | PR    | Note                                                                                                                                                                                       |
| ----- | ------------------------------------------------------------------ | ------------ | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| P1    | `ungzip` JS CRC after native inflate                               | fixed        | #2044 | native `DecompressionStream('gzip')` checks CRC/ISIZE ([zip.ts:216](../../packages/core/src/util/zip/zip.ts#L216)); 50 MB: 200 vs 1,520 ms (`gz.ts`, base A/B). New costs: B17                |
| P2    | `OrderedSet.getAt` re-dispatched in inner loops                    | still open   |       | [boundary.ts:44](../../packages/core/src/math/geometry/boundary.ts#L44) unchanged; `getAt2` 30 ms of `getBoundary` on flat-192                                                             |
| P3    | `subtractIS` binary search per element                             | fixed        | #2046 | one two-pointer pass ([ordered-set.ts:319](../../packages/core/src/data/int/impl/ordered-set.ts#L319)); tests pass                                                                         |
| P4    | `makeBuckets` `has`+`get`, second `getKey` pass                    | fixed        | #2045 | run reuse + ordinals ([buckets.ts:33](../../packages/core/src/data/util/buckets.ts#L33)); A/B too noisy at load 30 to re-measure; code checked. Small new cost: an `Int32Array(n)` even when already bucketed (4.5 MB transient on flat-192, per pass) |
| P5    | sorted models: `Table.view` copies every `atom_site` column        | fixed        | #2036 | lazy BCIF fields decode on `storage()`; dictionary string columns gather indices only                                                                                                      |
| P7    | `SortedRanges.firstIntersectionIndexFrom` linear from range 0      | fixed        | #2046 | binary search on range max ([sorted-ranges.ts:92](../../packages/core/src/data/int/sorted-ranges.ts#L92)); the per-residue cost left is in model `polymer.ts` (cross-slice note below)     |
| P8    | `isParamEqual` arrays via `Object.keys`                            | fixed        | #2046 | `shallowEqualArrays` / `arrayEqual` ([param-definition.ts:754](../../packages/core/src/util/param-definition.ts#L754))                                                                      |
| P9    | `distinctColors` uncached                                          | fixed        | #2046 | LRU(16) keyed by count + JSON of props; returns a copy. B6 is cached along with it                                                                                                         |
| P10   | `arrayView` gather polymorphic                                     | fixed        | #2032 | per-kind `pick*`; `aview.ts` (22 cols × 1.1M): `Column.view` 167–448 ms vs generic copy 716–2,858 ms in the same rounds                                                                   |
| P11   | NW: O(n·m) traceback unbounded in the model path; mixed matrix kinds | still open |       | [needleman-wunsch.ts:81](../../packages/core/src/algorithms/sequence/needleman-wunsch.ts#L81) `Uint8Array((n+1)·W)`; model `alignment.ts` still has no cell cap                              |
| P12   | `ChunkedArray.add*` stores polymorphic                             | still open   |       | unchanged ([chunked-array.ts:59](../../packages/core/src/data/util/chunked-array.ts#L59)); `chunk.ts` 3M `add3`: mono 28–48 ms, after 3 other ctors 35–90 ms (load 30)                     |
| B1    | `memoize1` never evicts; `getUniqueComponentNames` keyed by `BasicData` | still open |     | [memoize.ts:29](../../packages/core/src/util/memoize.ts#L29) a Map; [properties.ts:189](../../packages/model/src/formats/structure/basic/properties.ts#L189); plugin/ui `parameters.tsx:962/964` |
| B2    | typed column `toArray({ array })` ignores the ctor                 | partly fixed | #2036 | BCIF coordinates honor Float32; `arrayColumn` typed path ([column.ts:501](../../packages/core/src/data/db/column.ts#L501)) and masked BCIF `toFloatArray` still return the source kind (`colfuzz.ts`: the only mismatches) |
| B3    | `forEachSegment` last segment twice                                | fixed        | #2046 |                                                                                                                                                                                            |
| B4    | `LRUCache.set` duplicate key                                       | fixed        | #2046 | replaces; O(capacity) scan per `set` (fine at 16)                                                                                                                                         |
| B5    | `ChunkedArray.is` always falsy                                     | fixed        | #2046 |                                                                                                                                                                                            |
| B6    | `distinctColors` reads `samples[]` for `sampleList[]`              | still open   |       | [distinct.ts:157](../../packages/core/src/util/color/distinct.ts#L157), :162; now cached, same bias                                                                                       |
| B7    | `IntervalIterator`, `projectValue` wrong                           | fixed        | #2046 |                                                                                                                                                                                            |
| B8    | `getArrayBounds` clamps start to `rowCount - 1`                    | fixed        | #2046 | `colfuzz.ts`: 49k view/window/toArray-bound checks, 0 mismatches (besides B2)                                                                                                             |
| B9    | `deepClone` flattens typed arrays, Map, Set                        | fixed        | #2046 |                                                                                                                                                                                            |
| B10   | `ungzip` multi-member files throw                                  | still open   |       | both paths: `gz.ts` multi-member → "Trailing junk found after the end of the compressed stream" (also a `.gz` padded with trailing zeros, which `gunzip` accepts); base the same       |
| B11   | FEXTRA skip misses XLEN                                            | fixed        | #2032 |                                                                                                                                                                                            |
| B12   | failed asset read keeps its hold                                   | fixed        | #2033 | `_readHeld` releases on failure ([assets.ts](../../packages/core/src/util/assets.ts)); #2039 kept it                                                                                      |
| B13   | NW gap-tie doc                                                     | fixed        | #2046 | doc                                                                                                                                                                                        |
| B14   | `parseInt` wraps above 2³¹−1; `parseFloat` not correctly rounded   | still open   |       | [number-parser.ts:40](../../packages/core/src/util/number-parser.ts#L40); #2036 copied the same arithmetic into `parseIntBytes`/`parseFloatBytes` (:118), on purpose ("bit-identical")     |
| B15   | `_utf8Read` 2/4-byte fallback                                      | fixed        | #2046 |                                                                                                                                                                                            |
| B16   | `Column.view`/`window` of a masked column read raw slots           | fixed        | #2044 | io `binary/field.ts` storage is a masked copy; dictionary fields with a non-`''` masked row keep the string array                                                                         |

New bugs or costs from fixes since v2 (reviewed every hunk of the 57 changed files):

- #2044 (`ungzip` native path): allocates `ISIZE` from an untrusted trailer, now **twice** on failure (native, then the
  JS fallback), and a truncated file fails with an empty `TypeError` → **B17**. A corrupt-CRC file inflates twice
  (746 vs 200 ms on 50 MB) and still returns its data with a console warning (as v2).
- #2039 (`Sniffable` over bytes): `asciiIndexOf` rides `Uint8Array.prototype.indexOf`, a scalar loop in V8; the CIF
  variant guess on a plain mmCIF is now 2.2× the v2 string search → **P13**.
- #2043 (`value-palette`): a normalized palette without both ends builds a `Set` of every value as a string → **P14**.
- #2039 (`HeldAsset.file`): a new `File` copy of a download's bytes on every access → **P15**.
- #2045: `ordinals` allocated before the already-bucketed early return (transient 4 B/element; noted, not ranked).
- No correctness regressions found in `Column.storage`/dictionary views/windows (`colfuzz.ts`, 10 column kinds × 200
  random sizes), `subtractIS`, `firstIntersectionIndexFrom`, `IntervalIterator`, `projectValue`, `AsyncQueue`, LRU,
  `deepClone`, `utf8Read`/`utf8Encode`, `readResponseBody` (read; the core unit tests: 59 files, 543 tests pass).

## Open performance findings (ranked)

| Rank | Id  | Status     | Title                                                                                     | File                                                                                                   | Kind        | Estimate                                                                                     | Effort | Frozen |
| ---- | --- | ---------- | ----------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ | ----------- | -------------------------------------------------------------------------------------------- | ------ | ------ |
| 1    | P13 | new        | `asciiIndexOf` is a scalar `TypedArray.indexOf` scan; the CIF variant guess runs it 5× whole | [utf8.ts:271](../../packages/core/src/util/utf8.ts#L271), [provider.ts:190](../../packages/plugin/core/src/state/formats/provider.ts#L190) | time        | measured: 603 ms per `auto` open of a 124 MB mmCIF (268 ms in v2's string search; memchr would be 4 ms per scan) | S–M    | no     |
| 2    | P14 | new        | `numericRange` makes a string `Set` of every value to find min/max                        | [value-palette.ts:641](../../packages/core/src/util/color/value-palette.ts#L641)                       | time+memory | measured: 1.1M values 0.9–1.2 s, 90 MB alloc vs a min/max loop 37–96 ms                     | S      | no     |
| 3    | P12 | still open | `ChunkedArray.add`/`add3`/`add4` stores polymorphic across builders                       | [chunked-array.ts:59](../../packages/core/src/data/util/chunked-array.ts#L59)                          | time        | measured micro +25–90 % (load 30); v2's builder extrapolation stands (up to −100 ms flat-192 ball-and-stick) | M      | no     |
| 4    | P11 | still open | Needleman-Wunsch traceback unbounded in the model path; mixed matrix kinds                | [needleman-wunsch.ts:81](../../packages/core/src/algorithms/sequence/needleman-wunsch.ts#L81)          | time+memory | v2: 5k chain 25 MB / 135 ms, 20k 400 MB / 2.1 s; +27 % mixed kinds (code unchanged)          | M      | no     |
| 5    | P2  | still open | `OrderedSet.getAt` re-dispatch per element in `getBoundary`                               | [ordered-set.ts:37](../../packages/core/src/data/int/impl/ordered-set.ts#L37)                          | time        | v2: up to 25 % of small-entry boundary time (7qpd 1.7 ms); flat-192 30 ms of 523 ms          | M      | no     |
| 6    | P15 | new        | `HeldAsset.file` copies a download's bytes into a new `File` on every access              | [assets.ts:155](../../packages/core/src/util/assets.ts#L155)                                           | memory      | estimated: one extra transient copy per held download per session save (`snapshots.ts:437`) | S      | no     |

### P13: byte search is a scalar loop, five times per guess

- Mechanism: `Sniffable.searchBytes` → `asciiIndexOf` finds candidates with `Uint8Array.prototype.indexOf(anchor)`. V8
  implements typed-array `indexOf` as a per-element loop (no memchr); Node's `Buffer.indexOf` is memchr.
  `sniffCifVariant` asks five needles of every `.cif` in `auto` (`_volume_data_3d_info.`, `_refln.pdbx_FWT`,
  `_refln.pdbx_DELFWT`, `atom_site_fract_x`, `atom_site.fract_x`); a plain mmCIF has none, so each scans all of it
  (anchor candidates are few: 7–451, so the time is the scan). The `lastSniff` cache makes it once per data.
- Evidence (`sniff4.ts`, 124 MB flat-192.cif, min of 7, load ~27): `asciiIndexOf` × 5 603 ms; `string.indexOf` × 5
  on the same text 268 ms; one `Uint8Array.indexOf(absent byte)` 95 ms; `Buffer.indexOf` 4 ms. v2 searched the
  decoded text (`includes`), so the `auto` open of a large text CIF pays about +0.3 s now, before parsing.
- Fix: one pass for all needles of the guess (they name 3 categories), or search category names only where they can
  be (a line starting with `_`); for the scan itself, a 32-bit SWAR word scan for the anchor (`sniff3.ts`: 1.1–4.0 s
  vs 1.9–6.1 s for `asciiIndexOf`, same rounds) or decoding the first and last MB only if the formats allow it. Risk:
  low; keep `asciiIndexOf` semantics (byte index).

### P14: normalized domain through strings

- Mechanism: `numericRange` → `distinct(values(), false, cifMissing)` converts every value to a string into a `Set`,
  spreads it, `parseFloat`s each back, then two `reduce`s. Callers: Seq* `makePaletteFunction` passes every row value
  (`palette.ts:125`); MVS passes distinct values already (MVS itself not re-audited here).
- Evidence (`vpal.ts`, 1.1M B-factor-like values, ~33k distinct): 909–1,219 ms and 90–95 MB allocated per palette;
  a plain min/max loop over the same values 37–96 ms.
- Fix: a loop with `toNumber(value, cifMissing)` skipping NaN and missing (same result: `parseFloat(\`${n}\`)` is `n`
  for numbers; booleans are NaN both ways). S.

### P15: `file` getter copies

- `get file()` returns `this.source ?? new File([this.bytes!], 'raw-data')`: every access copies the bytes into a Blob.
  Session save (`plugin/core/state/manager/snapshots.ts:437`) reads `file` then `arrayBuffer()`, two transient copies
  per held download where v2 had one. Fix: cache the made `File` until the bytes change, or give the serializer the
  bytes. Estimated, not measured.

### P12, P11, P2

Unchanged code; v2's evidence and fixes stand ([../perf-audit-2/01-core-data-util.md](../perf-audit-2/01-core-data-util.md)).
P12 re-run (`chunk.ts`) still shows the polymorphic slowdown; P2 and P11 not re-measured (code identical).

### Cross-slice note (slice 05): always-true intersection per residue

[polymer.ts:25](../../packages/model/src/model/structure/structure/util/polymer.ts#L25) tests
`OrderedSet.areIntersecting(Interval.ofRange(elements[start], elements[end - 1]), elements)` for every residue
segment of `elements`; `elements[start]` is in `elements`, so it is always true. flat-192 structure profile:
`areIntersecting*` 37 ms + `Interval.ofRange` 21 ms of `getAtomicPolymerElements`' 94 ms (load ~30). Fix: drop the test.

## Bugs

| Id  | Severity | Status     | Title                                                                                      | File                                                                                       | Failure scenario                                                                                                   |
| --- | -------- | ---------- | ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------ |
| B17 | medium   | new        | `ungzip` sizes its buffer from the trailer: a truncated `.gz` allocates up to 4 GiB, twice  | [zip.ts:216](../../packages/core/src/util/zip/zip.ts#L216), [zip.ts:266](../../packages/core/src/util/zip/zip.ts#L266) | a 5 MB `.gz` cut short → 1.0–3.4 GiB `arrayBuffers`, RSS +2.4–6.9 GB in some runs; then an empty `TypeError`       |
| B1  | medium   | still open | `memoize1` never evicts; `BasicData`-keyed caller pins each load                            | [properties.ts:189](../../packages/model/src/formats/structure/basic/properties.ts#L189)   | every mmCIF load without `chem_comp` keeps its source (v1: +2.7–4.1 MB per load)                                   |
| B10 | low      | still open | multi-member (and zero-padded) gzip throws                                                 | [zip.ts:210](../../packages/core/src/util/zip/zip.ts#L210)                                 | `gzip -c a b > ab.gz` → "Trailing junk found after the end of the compressed stream"                               |
| B14 | low      | still open | `parseInt` wraps above 2³¹−1; `parseFloat` not correctly rounded (now also the `*Bytes` twins) | [number-parser.ts:40](../../packages/core/src/util/number-parser.ts#L40)                  | `'3000000000'` → −1294967296; `'1.805'` → 1.8050000000000002                                                       |
| B2  | low      | partly     | typed `toArray({ array })` ignores the ctor                                                 | [column.ts:501](../../packages/core/src/data/db/column.ts#L501)                            | `Column.view(f32Col, map).toArray({ array: Float64Array })` returns Float32Array (`colfuzz.ts`)                    |
| B6  | low      | still open | `distinctColors` clustering reads `samples[]` for `sampleList[]`                            | [distinct.ts:157](../../packages/core/src/util/color/distinct.ts#L157)                     | 'generate' palettes biased; now cached                                                                             |
| B18 | low      | new        | `Column.areEqual` of two array columns compares storage before the schema's transform      | [column.ts:635](../../packages/core/src/data/db/column.ts#L635)                            | a lowercase str column over `['A']` vs one over `['a']`: values equal, `areEqual` false; dictionary columns compare values (inconsistent) |

- **B17**: `ungzip` (native path, #2044) does `new Uint8Array(gzipDeclaredSize(file))`; on failure `ungzipBlocks`
  allocates `new Uint8Array(isize)` again. In a truncated file the "ISIZE" is 4 bytes of deflate data: up to 4 GiB.
  `trunc2.ts` (5 MB of text gzipped, cut by 1–5,000 bytes): declared 362–3,438 MiB; `process.memoryUsage().arrayBuffers`
  372–3,448 MiB after the call; RSS +2.4, +3.3, +6.9 GB in 3 of 9 cuts (HEAD), v2 base the same order on one
  allocation. Error: `TypeError` with an empty message (`trunc.ts`), the same in v2. `assets.ts` guards with
  `MaxDeflateRatio` (:125) only for `decompress` reads; `data-source.ts:164` (a `.gz` file opened by name),
  `plugin/core/.../data/fetch.ts:240` and `prism/viewer/.../host.ts:504` call `ungzip` unguarded. In a browser a
  multi-GiB `ArrayBuffer` either throws `RangeError` (then the fallback throws the same) or commits memory on a tab
  with a ~4 GiB budget. Fix: move the ratio guard into `ungzip`; above it (or always on the native path) inflate into
  chunks without a pre-sized buffer; skip the JS fallback when the native error is not "multi-member"-like; wrap the
  error with "corrupt or truncated gzip".
- **B18**: pre-existing semantics, now inconsistent with the new dictionary branch (:637). Only caller:
  `Table.areEqual` in `atomic.ts:158` (sharing the hierarchy between frames), so a false negative costs memory, not
  correctness. Fix: compare storage only when the schema has no transform, else values.
- B1, B10, B14, B2, B6: v2's evidence and fixes stand.

## Already good, re-checked

- `Column.storage` + dictionary columns: views, view-of-view, windows (empty, full, partial), `toArray` bounds and
  `Array` ctor, mixed-type string arrays, transforms: 49,234 checks, 0 mismatches (`colfuzz.ts`).
- `fetchData`/`readResponseBody`: a known `Content-Length` reads into one buffer; overflow and bogus lengths fall back
  to chunks; abort listener removed in `finally`.
- `utf8Read` short-ASCII path, `utf8ReadValue` BOM, `utf8Encode` surrogate split at piece ends.
- `AsyncQueue` per-value waiters; `LRUCache.set` replace; `deepClone` typed arrays/Map/Set.
- `text-measure.ts` (bounded cache, cleared past 65,536 widths); `value-palette` `at` paths allocation-free.

## Not audited

- `value-palette-params.ts`, `split-color.ts` parsing in depth; MVS callers (covered by `.v6/mvs-optimizations/`).
- Wall-clock re-measurement of P4 (`buckets.ts`), P2, P11: load 26–33 made A/B differences unreadable; status by code.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/01/`;
build with `node build.mjs $PWD/x.ts $PWD/x.mjs` (HEAD) or `node build-base.mjs ...` (v2 base in `./base`).

| Script                    | Settles                                                         |
| ------------------------- | --------------------------------------------------------------- |
| `aview.ts`                | P10 fixed: `Column.view` vs generic and per-kind copies         |
| `chunk.ts`                | P12 still polymorphic                                           |
| `gz.ts`                   | P1 fixed (base A/B), B10 multi-member and zero padding, bad CRC |
| `trunc.ts`, `trunc2.ts`   | B17 empty error, declared sizes, `arrayBuffers`, RSS            |
| `colfuzz.ts`              | Column views/windows/toArray, B2, B18 (`B2=1`, `EQ=1`)          |
| `sniff3.ts`, `sniff4.ts`  | P13 `asciiIndexOf` vs string search, memchr, SWAR, one-pass     |
| `vpal.ts`                 | P14 normalized domain                                           |
| `buckets.ts`              | P4 A/B (inconclusive under load)                                |
| `prof192/`                | flat-192 parse, model, structure, lookup+boundary profiles      |
