# CIF parsing (text and BinaryCIF)

How `@molstar/io` reads CIF, and why. `packages/io` has no `_spec`; this note is its record for the CIF readers. Numbers:
Apple M4, Node 25, `scripts/pipeline-bench.mjs` (flat-192 = 1tup ×192 as one model, 1.12M atoms).

## Text CIF: a byte tokenizer

- **The tokenizer reads bytes** (`reader/cif/text/parser.ts`): a `Uint8Array` of UTF-8. Delimiters (whitespace, quotes,
  `;`, `_`, `#`) are ASCII and never occur inside a UTF-8 multi-byte sequence, so scanning bytes is safe; values decode
  as UTF-8 when read. Why: `charCodeAt` over a string is ~2× slower than indexing bytes, and slower still when the call
  site has seen a `ChunkedBigString` (flat-192-cif parse 594–726 → 285–327 ms).
- **A string is encoded up front** (`utf8Encode`, core): one `encodeInto` pass for ASCII (~15 ms for 124 MB), then the
  one tokenizer. Not a second tokenizer over strings: two would diverge, and the copy costs less than the string scan
  saves. The cost is memory while the caller keeps its string: the parsed file holds the bytes (124 MB on flat-192) where
  it held a reference to the string.
- **Tokens** (`text/field.ts`, internal): per loop, one column-major pair of arrays, a Uint32 start and a Uint8 length
  per value (255 = look the end up in a map, for values of 255+ bytes). 5 B a value against 8 for a start and an end;
  one pair per loop keeps the reading loop on locals, column-major keeps a field's values contiguous for its readers. A
  loop grows from the density read so far and is trimmed. flat-192: 132 MB of tokens against 320 MB before (8 B plus
  the 1.61× growth slack). The shared `Tokens`/`Tokenizer` of the other text formats are unchanged.
- **`str`** decodes per value: ASCII up to 24 bytes with `String.fromCharCode` (`utf8ReadValue`, core; a `TextDecoder`
  call is ~0.4 µs; a leading U+FEFF of a value is kept), and a field of 64 rows or more has a 256-slot cache of values
  up to 8 bytes (up to 4 keyed by the bytes), so repeated names are made once and keep their computed hash. Fewer rows
  have none: its ~3 KB outweighed a dictionary's single items (the CCD). A category of single items shares one token
  store. Numbers parse on bytes (`parseIntBytes`, `parseFloatBytes`, core), bit-identical to the
  string parsers.
- **Compatibility**: every field accessor equals the string tokenizer's (kept as the test oracle,
  `text/_test/string-parser.ts`); error lines are the same; a UTF-8 byte order mark is skipped (a string with one failed
  before). Input of 4 GiB or more is refused (Uint32 starts). A lone surrogate in a string input becomes U+FFFD (the
  encode). A file of single items only (100k blocks of 7) parses ~18 % slower than with strings, a typed-array pair per
  category; with loops (the CCD's atoms and bonds) it is faster.

## `CIF.parse` and the plugin

- `parseCifText(data: StringLike | Uint8Array)`: additive.
- `CIF.parse(Uint8Array)` meant BinaryCIF; it now looks at the first byte: a MessagePack map marker (`0x80–0x8f`,
  `0xde`, `0xdf`) is BinaryCIF, anything else text (a UTF-8 byte order mark is skipped). A BinaryCIF file always starts
  with a map marker, so only input that failed before (text CIF as bytes) reads differently. The plugin's `ParseCif`
  goes through it, so a Binary data cell holding text CIF parses.
- Bytes that are neither (gzip, say) now fail as text, "Unexpected token", where they failed as MessagePack; whitespace
  only parses to no blocks. No bytes fail as BinaryCIF, as before.
- A file parsed from bytes reads them in place: a caller that reuses the buffer changes the values (a string could not
  change).
- `CIF.parseBinary` stays strict, no sniffing: `ParseBlob` calls it for its binary entries.
- The plugin keeps reading `.cif` and `.mmcif` files as strings (the maintainers' decision): reading them as binary
  would skip the encode copy and the 512 MB string limit, but would change the Data node's type of opened files. The
  way to get both later is a byte-backed `StringLike` (`perf-audit-2/13-lammps-large-text.md`).

## BinaryCIF

- **A field decodes on first access** (`binary/field.ts`), given its category's row count; `createColumn` forwards
  `storage()` on call, and a sorted or windowed `atom_site` makes each column on its first read (`model` `basic/sort.ts`),
  so the columns a model never reads are neither decoded nor copied (perf-audit-2 04/P5).
- **Fused decoding** (`common/binary-cif/decoder.ts`): `[FixedPoint | IntervalQuantization] > [Delta] > [RunLength] >
  [IntegerPacking] > ByteArray` with two or more steps (the chains a survey of PDBe/RCSB files found: coordinates, B,
  occupancy, ids, sequence ids, string indices and offsets) decode in one pass with no intermediate arrays; a
  StringArray's index chain maps straight to its strings. Bit-identical to the step decoder (`decodeSteps`, the oracle
  of a randomized test): values wrap as the stores into the Int32 intermediates do. Other chains take the steps. Two
  loop copies, numbers and strings: one function writing typed and JS arrays went megamorphic.
- **Float32 straight**: a quantized column read with `toFloatArray({ array: Float32Array })` decodes to Float32 (each
  value computed in double, rounded once); the other accessors decode to the source type, so writers see exact
  values. Sorted and windowed models gather the coordinates from that array, of its type, reading only their own rows. The
  conformation of a quantized BinaryCIF model (the PDB's coordinates) is therefore Float32, like a text one's; it was
  Float64 because `typedArrayWindow` ignores `array`. A raw Float64 column stays Float64 in every path. A window of a
  typed column asks its column for another array kind, so IHM coarse coordinates (`Table.window` per model) are Float32
  for one model and many alike.
- **String columns as dictionaries** (04/P16): a StringArray column decodes to its distinct strings and an index per
  row in the smallest unsigned type (`decodeStringDictionary`), not a string per row; its storage is
  `{ kind: 'dictionary', indices, values }` (`Column.storage`, core). Views and windows gather or slice the indices,
  `Column.areEqual` compares them, and `sortAtomSite` buckets by them when the strings are distinct. The values are
  those of the array of strings: an index out of range, or a masked row whose string is not `''`, keeps the array.
  flat-192 model stage: 419 → ~334 ms, retained 244 → 124 MB.
