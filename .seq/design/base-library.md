# The base library: reuse, extend, promote, add

## Intent

The author, 2026-11-17: "in general, should we get a catalog of precetents to improve on our base library?"

The root `AGENTS.md` says: reuse before writing helpers; a generic routine missing from `@molstar/core` is proposed
there with a test, not kept as a private copy. The previous implementation kept many private copies; this catalog says,
for each need of Seq\*, what the base library already has, what to extend, what to promote from the previous code, and
what is missing. External precedents are [sources/precedents](../sources/precedents/README.md).

## Provisional

Every "extend", "promote" and "add" below changes a public API of `@molstar/core`, `@molstar/io` or `@molstar/ui`, so
each is the author's call when its work starts. "Reuse" needs no decision. Promoting means taking the previous code as
a guide and writing the generic piece in core with its test; `packages/seq` never imports the previous packages.

## Reuse as is

| Need | Where | Names |
|---|---|---|
| hashing, stable hashes of values | `@molstar/core/data/util/hash-functions` | `hashFnv32a`, `murmurHash3_128`, `hashMurmur128o`, `hashFnv256o` |
| integer sets and ranges | `@molstar/core/data/int/*` | `Interval`, `SortedArray`, `OrderedSet`, `SortedRanges`, `Segmentation` |
| sorting | `@molstar/core/data/util/sort` | `sortArray`, `radixSort`, `radixSortIndices` |
| chunked and scratch arrays | `@molstar/core/data/util/chunked-array`, `@molstar/core/util/int-buffer` | `ChunkedArray`, `IntBuffer`, `ScratchMemory` |
| memoization, reference caches | `@molstar/core/util/memoize`, `reference-cache` | `memoizeLatest`, `memoize1`, `createReferenceCache` |
| assets: downloads and files, ref-counted, one download shared, gzip inflated once | `@molstar/core/util/assets` | `AssetManager`, `Asset.Url`, `Asset.File`, `resolve(asset, type, { shared, decompress })` |
| tasks, cancellation, progress | `@molstar/core/task` | `Task`, `RuntimeContext`, `Progress`, `chunkedSubtask`, `Scheduler` |
| zip, gzip, inflate, deflate | `@molstar/core/util/zip/zip` | `unzip`, `zip`, `ungzip`, `inflate`, `deflate` (limits included) |
| MessagePack | `@molstar/core/data/msgpack` (moved down from io, which re-exports it; no allowlist entry needed) | `encodeMsgPack`, `decodeMsgPack` with their extension hooks: the repr codec with typed columns (`@molstar/core/repr/codecs/msgpack`) is built on them ([tree-core.md](tree-core.md)) |
| global alignment | `@molstar/core/algorithms/sequence/needleman-wunsch` | `globalAlign`, `GlobalAlignScoring` |
| color, contrast, lists, scales | `@molstar/core/util/color/*` | `Color.luminance`, `Color.contrast`, `ColorLists` (Okabe–Ito, Cividis), `ColorScale`, `distinctColors` |
| histograms, array statistics | `@molstar/core/math/histogram`, `@molstar/core/util/array` | `calculateHistogram`, `arrayMinMax`, `arrayStats` |
| number and value formatting | `@molstar/core/util/number`, `string-format` | `toPrecision`, `FormatTemplate` |
| params, bindings | `@molstar/core/util/param-definition`, `binding` | `ParamDefinition`, `Binding` |
| persistent collections | `@molstar/core/data/immutable` (new, [immutable.md](immutable.md)) | `HashMap`, `HashSet`, `OrderedHashMap`, `OrderedHashSet`, `Vector`, transients, cursors |
| reactive state | `@molstar/core/reactive` | `State`, `derive`, `KeyedState`, `ReactiveModel` |
| history, commands, keymaps | `@molstar/ui/history`, `commands` | `HistoryPort` (checkpoints) as the adapter behind Seq\*'s own history port, in a composition root only ([changes.md](changes.md)); `Keymap`, `CommandManager`; `forms` later, with the settings form |

## Extend in place

| Need | Where | Extension |
|---|---|---|
| text layout | `@molstar/core/util/text-measure` | pretext's prepare and layout split, lines for the text view, inline boxes (chips); measure function as the port's adapter ([rendering.md](rendering.md)) |
| byte-bounded caches | `@molstar/core/util/lru-cache` | O(1) lookup, a bound in bytes beside the count (the source cache, the glyph atlas) |
| fetching by range | `@molstar/core/util/assets` (`AssetManager`; its `ResolveOptions` carry credentials, redirects, the signal and a streamed byte limit, and `checkAssetManager` is the conformance suite every implementation passes: security.md) | range reads of an asset: a page cache by asset and page, ref-counted like its entries, shared across readers, aborted by the task; `Blob.slice` for files (the asset port's default adapter, [ports.md](ports.md)) |
| BGZF and tabix | `@molstar/io/common/bgzf`, `common/tabix` (the fifth slice, [genome-formats](genome-formats.md), beside io's readers; no `.gzi`); core's `ungzip` over a multi-member stream | blocks and virtual offsets, TBI and CSI indexes, a region's lines and `regionBytes`; every bound refused before it allocates |
| one luminance | `@molstar/ui/stage/canvas-ink` | use `Color.luminance` instead of its own |

## Promote: generic pieces the previous implementation kept private

Written anew in core with their tests, the previous code read as a guide.

| Piece | Goes to | What |
|---|---|---|
| half-open range sets | `@molstar/core/data/int` | sorted, merged half-open ranges in `Int32Array`s with union, intersection, subtraction and overlap by binary search, taken from Seq\*'s Location (the author, 2026-11-17: "ok") |
| interval index | `@molstar/core/data` | an implicit interval tree over sorted typed columns (overlap queries by region) |
| summary pyramids | `@molstar/core/math` | value pyramids (min, max, mean, count) and class pyramids by level ([performance.md](performance.md)) |
| genetic codes | `@molstar/core/algorithms/sequence/` (a module beside `needleman-wunsch`) | NCBI codon tables, Sec and Pyl, one- and three-letter codes |
| numberings | `@molstar/core` (beside the tree core, or in sequence algorithms) | offset, segments, labels with insertion codes, names; label order |
| residue alphabets and schemes | `@molstar/core/util/color` and sequence constants | residue color schemes as rules (Clustal X, Zappo, Taylor, ESPript), alphabets beside `@molstar/model`'s one-letter codes (moved down to core) |
| MessagePack | `@molstar/core/data/msgpack` (the encoder and decoder) and `@molstar/core/repr/codecs/msgpack` (the tree's codec), moved down from `@molstar/io`, which re-exports both | JSON's binary twin for sources, read by `seqx-core` without an io dependency (the author, 2026-11-19; built with FS) |
| a delimited-record tokenizer | `@molstar/core/data/text/delimited` | RFC 4180 records with plain string cells (a one-character delimiter, quote and comment, a BOM skipped), written from io's CSV parser as the guide; io's reader unchanged ([files.md](files.md), agreed 2026-11-19) |
| standard components | `@molstar/core/util/sequence-constants` (built in the fourth slice, ML1: `SequenceAlphabet`, `StandardComponents`, `standardComponent`, `standardLetter`, and `AminoAcidNames` moved there as the one table) | the standard component of each letter by alphabet (A ↔ ALA, DA, A), the inverse half of the row above ([molecules](molecules.md), agreed 2026-11-19) |
| substitution matrices | `@molstar/core/algorithms/sequence/` | BLOSUM and friends, now in `@molstar/model` |

## Add: missing everywhere

| Piece | Goes to | Why |
|---|---|---|
| the tree core | `@molstar/core/repr` | decided ([tree-core.md](tree-core.md)) |
| an S-expression reader and printer | `@molstar/core/repr/codecs/sexpr` (with the tree core) | meso\*'s, generic |
| a worker pool | `@molstar/core/task` | when the compute port comes ([ports.md](ports.md)) |
| a ring buffer | `@molstar/core/data/util` | input and frame history without allocation |
| genomic unit formatting and tick steps | `@molstar/core/util/number` | kb and Mb, nice tick steps for rulers |
| a fixed seeded generator | `@molstar/core/data/util/hash-functions` (`PCG`) | decided, the author, 2026-11-19: "Fix PCG in core": today's `PCG` keeps 32 bits of state with a 64-bit output step and is biased (CS: G 33%, C 20% over 1 Mb; consecutive draws correlated). A real PCG32 (64-bit state, its `int()`/`float()` kept); Mol*'s `distinct.ts` colors and every seeded Seq\* fixture regenerated once, reviewed as a diff; a behavior change in `breaking-v6-changes.md`. Built (PG): its tests match the reference's outputs and a BigInt oracle of the published algorithm over many seeds, the rare low-to-high carry included (§seeded-generators-checked) |
| SIFTS tables | `@molstar/io` | SIFTS residue tables read outside mmCIF, for the structure map's source ([maps.md](maps.md)) |

## Duplicates to remove on the way

- Two luminance functions (`Color.luminance` and `canvas-ink`): one.
- Two schedulers (`@molstar/core/task/util/scheduler` and the frame scheduler the viewer needs): the viewer's frames go
  through `@molstar/ui/dom`'s `DomHost`, not a third scheduler.
- The previous implementation's own sets, colors, hashes, inflate and text measuring: not carried over; Seq\* uses
  core's.

## Consequences

- `packages/seq` has no private copy of a generic routine; a piece written here first and found generic moves to core
  with its test.
- Each extension, promotion or addition is its own small change to the base package, with its test, its `_spec` line
  and the v6 changelog entry, done when the Seq\* work first needs it.
