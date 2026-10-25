# CIF

Text CIF (mmCIF, CCD, BIRD, dictionaries, STAR) and BinaryCIF: `src/reader/cif.ts` (`CIF`), `src/reader/cif/`,
`src/common/binary-cif/`, `src/writer/cif*`. Both forms parse to the same data model; readers above it (model's mmCIF
and CIF-core readers, volume server data, MVS annotations) never see which.

## 1. Data model (`reader/cif/data-model.ts`)

| Type          | What                                                                                                                         |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `CifFile`     | `blocks`, an optional `name`                                                                                                 |
| `CifBlock`    | a `CifFrame` (`header`, `categoryNames` in file order, `categories` by name) plus `saveFrames`; `getField('cat.field')`      |
| `CifCategory` | `name` (no leading `_`), `rowCount`, `fieldNames`, `getField(name)`                                                          |
| `CifField`    | one column: `str`/`int`/`float`/`valueKind` per row, `areValuesEqual`, `to*Array`, `isDefined`, `binaryEncoding`, `storage?` |

- **Names are case-insensitive for lookup**: `getField` tries the exact name, then the lowercase one (text and binary
  categories both keep a lowercase map). Category keys are as written.
- **`.` and `?`**: `valueKind` is `NotPresent` / `Unknown`; `str` gives `''`, `int`/`float` 0. Text CIF does not tell
  the bare tokens from the quoted values `'.'` and `'?'` (the parser's header lists its other deviations from CIF 1.1:
  keywords only are case-insensitive, no `\` line folding).
- **A field's functions do not use `this`**: callers pass `field.str` around unbound (`schema.ts`, model readers).
- `CifField.ofStrings/ofNumbers/ofColumn/ofTokens/ofTextByteTokens`, `CifCategory.ofFields/ofTable`: data made in code
  (writers, converters such as PDB → CIF) in the same shape.

## 2. `CIF.parse` (`reader/cif.ts`)

- **A string is text CIF; bytes are told by their first byte**: a MessagePack map marker (`0x80–0x8f`, `0xde`, `0xdf`)
  is BinaryCIF (a BinaryCIF file is one map), anything else text (`data_`, `#`, whitespace, a byte order mark). Empty
  bytes go to the binary parser and fail there. Why: the plugin holds downloads and files as bytes whatever they are,
  and one entry point serves `ParseCif`.
- **`CIF.parseText`** (`parseCifText`) takes a `StringLike` or bytes; **`CIF.parseBinary`** is strict BinaryCIF, no
  sniffing (a `Data.Blob`'s binary entries).
- Bytes that are neither (gzip) fail as text ("Unexpected token"); whitespace only parses to a file of no blocks.
- `CIF.schema.*` and `toDatabase` give typed views (§5).

## 3. Text CIF (`reader/cif/text/`)

- **The tokenizer reads bytes** (`parser.ts`); a string is encoded first (`utf8Encode`). A leading UTF-8 byte order mark
  is skipped. Input of 4 GiB or more is refused (Uint32 starts). A lone surrogate in a string input becomes U+FFFD (the
  encode). The reasons are the package's (README, Decisions).
- **Tokens** (`field.ts`, internal, `ByteTokens`): per loop one column-major pair of arrays, a Uint32 start and a Uint8
  length per value; length 255 means 255 bytes or more, the end is in a map (`long`). 5 B a value against 8 for start
  and end; one pair per loop keeps the reading loop on locals; column-major keeps a field's values contiguous. A loop's
  arrays grow from the density read so far and are trimmed at its end. A category of single items (no `loop_`) shares
  one token store for all its fields. Loops are read in chunks (`chunkedSubtask`): progress and abort between them.
- **`str`** decodes per value (`utf8ReadValue`, core): ASCII up to 24 bytes with `String.fromCharCode` (a `TextDecoder`
  call costs ~0.4 µs); a leading U+FEFF inside a value is kept. A field of 64 rows or more reads values of up to 8 bytes
  through a `ByteStringCache` (256 slots; up to 4 bytes keyed by the bytes as an int, 5–8 by start, compared by bytes):
  repeated names (atoms, residues, chains, elements) are one string each, with V8's hash computed once. The cache turns
  itself off for a column whose first 4096 lookups mostly missed (unique ids pay only the miss until then). Under 64
  rows there is no cache: its ~3 KB outweighs a dictionary's single items (the CCD).
- **Numbers** parse from the bytes (`parseIntBytes`, `parseFloatBytes`), bit-identical to the string parsers.
- **The oracle**: every field accessor equals the string tokenizer's (`text/_test/string-parser.ts`, kept for tests
  only), and so do the error lines. A file of single items only (100k blocks of 7) parses a little slower than over a
  string (a typed-array pair per category); files with loops are faster.
- **The parsed file keeps the bytes** and reads values from them: a caller that reuses the buffer changes the values
  (README, contract 3). A string input is held as its encoded bytes.

## 4. BinaryCIF (`reader/cif/binary/`, `common/binary-cif/`)

- **Format** (`encoding.ts`, version 0.3): MessagePack of blocks → categories → columns, each column `data` and an
  optional `mask`, each an `EncodedData` (bytes plus the encoding chain to undo: `ByteArray`, `FixedPoint`,
  `IntervalQuantization`, `RunLength`, `Delta`, `IntegerPacking`, `StringArray`). `parseCifBinary` checks the version
  (major and minor each at least 0.3) and decodes the MessagePack (`common/msgpack/decode.ts`; `bin` values are
  exact-size copies, so a decoder may view `data.buffer`).
- **A field decodes on first access** (`binary/field.ts`), given its category's row count; `getField` caches the field,
  `createColumn` (§5) forwards `storage()` on call. A model reads only the columns it uses; a sorted or windowed
  `atom_site` makes each column on its first read (model `formats/structure/basic/sort.ts`).
- **Fused decoding** (`decoder.ts`, `decodeFused`): the chains
  `[FixedPoint | IntervalQuantization] > [Delta] > [RunLength] > [IntegerPacking] > ByteArray` with two or more steps
  (coordinates, B, occupancy, ids, sequence ids, string indices and offsets in PDBe and RCSB files) decode in one pass
  with no intermediate arrays; a StringArray's index chain maps straight to its strings. Bit-identical to the step
  decoder (`decodeSteps`, the oracle of a randomized test): values wrap as the stores into the Int32 intermediates do.
  Other chains decode step by step. Two loop copies, numbers and strings: one function writing typed and JS arrays went
  megamorphic.
- **Float32 straight**: a quantized column read with `toFloatArray({ array: Float32Array })` decodes to Float32 (each
  value computed in double, rounded once: equal to a Float32 copy of the source-type array). Other accessors decode to
  the source type, so writers see exact values. Sorted and windowed models gather coordinates from that array, of its
  type, reading only their rows, so a quantized BinaryCIF model's conformation is Float32 like a text one's. A raw
  Float64 column stays Float64. A window of a typed column asks its column for the requested array kind (IHM coarse
  coordinates, `Table.window` per model).
- **String columns are dictionaries** (`decodeStringDictionary`): a StringArray column decodes to its distinct strings
  and one index per row in the smallest unsigned type; its storage is `{ kind: 'dictionary', indices, values }`
  (`Column.storage`, core). Views and windows gather or slice the indices, `Column.areEqual` compares them, and model's
  `sortAtomSite` buckets by them when the strings are distinct. The values are those of the array of strings: an index
  out of range, or a masked row whose string is not `''`, keeps the array instead.
- **Header without data** (`header.ts`, `getBinaryCifHeader`, memoized per buffer): the encoder, blocks, categories and
  column names, walked over the MessagePack without decoding columns: what format providers sniff by
  (`binaryCifHasCategory`, `binaryCifHasColumn`).

## 5. Schemas (`reader/cif/schema.ts`, `schema/*`)

- `toDatabase(schema, frame, aliases?)` and `toTable(schema, category)` wrap categories as typed `Table`s of `Column`s
  over the fields: no copy, values read through the field. A table's column is made on its first property access and
  cached (`CategoryTable`), so the columns a reader never touches cost nothing. A missing category is an empty table
  (`_isDefined` false), a missing field `Column.Undefined` (`isDefined` false).
- Schemas: `mmcif.ts` (code-generated by molstar/ciftools from the mmCIF, IHM and MA dictionaries; do not edit by hand;
  `mmcif-extras.ts` for the rest), `ccd.ts`, `bird.ts`, `dic.ts` and `ddl.ts`, `cif-core.ts` (with `CifCore_Aliases`: a
  field read by any of its dictionary names), `density-server.ts`, `segmentation.ts`, `sf.ts`.
- String columns apply the schema's `transform` (`lowercase`, `uppercase`) to `value` and `toArray`; `storage` stays the
  raw field's, so consumers of storage see untransformed values. List columns split `str`; tensor columns read the
  fields `key[i]…` or `key_ij…` (either naming, 0- or 1-based).

## 6. Writing (`writer/cif.ts`, `writer/cif/`)

- `CifWriter.createEncoder({ binary })`: a `TextEncoder` (string) or `BinaryEncoder` (bytes). A writer declares
  categories as `Field`s over a data source (`CifWriter.fields()`, `categoryInstance`); the encoder walks them, so data
  is written without being made into tables first.
- Binary encodings per field: given (`Field` options, `CifWriter.Encodings`: delta-RLE, fixed point), from a provider
  (`createEncodingProviderFromCifFrame`: reuse a parsed file's), or classified per array (`binaryAutoClassifyEncoding`,
  `classifier.ts`). `ArrayEncoder` (`array-encoder.ts`) implements the steps the decoder undoes.

## 7. Callers outside io

- **Format sniffing** reads the bytes without decoding them: `guessCifVariant` (plugin `state/formats/provider.ts`) and
  `looksLikeMmcifParticles` (model `formats/particles/mmcif-sniff.ts`), through core's `Sniffable` and `asciiIndexOf`
  (anchored on a needle's first lowercase letter, rare in CIF data rows). `guessCifVariant` keeps its last answer for
  the same bytes (every CIF provider asks in turn) and searches needles that share a prefix once, for the prefix: a
  plain mmCIF holds none of them, so every search reads the whole file.
- **The plugin parses what it holds**: `Download`, `ReadFile` and `DownloadBlob` give Response Data for `.cif` and
  `.mmcif`, and `ParseCif` (and `ParseBlob` for a response blob) hands its bytes to `CIF.parse`, shared with the asset
  manager (core `util/_spec/assets.md` §4): no decode to a string and encode back, one copy of the bytes, no 512 MB
  string limit. A byte-backed `StringLike` was the alternative, kept off: a second `StringLike` kind at every text
  parser's call sites.
