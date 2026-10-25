# Mol\* io: the specification

`@molstar/io` (`packages/io/`) parses and writes file formats as plain data: a parsed file (a `CifFile`, a format's own
record of tables and columns), never a model. Building models, volumes and shapes from it is `@molstar/model`'s
(`formats/`). It depends on `@molstar/core` only. This spec is for working on the code without reading all of it: where
things are, the contracts a reader keeps, and the decisions (the rule, then why). Code comments cite sections
(`cif.md §3`). The code is the truth for how.

| File             | What it covers                                                                             |
| ---------------- | ------------------------------------------------------------------------------------------ |
| [cif.md](cif.md) | CIF: the data model, `CIF.parse`, the text tokenizer, BinaryCIF decoding, schemas, writing |

## Map

| Path                                            | What                                                                                                               |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `src/reader/cif.ts`, `src/reader/cif/`          | CIF and BinaryCIF readers, `data-model.ts` (`CifFile` … `CifField`), `schema.ts` and `schema/*` (cif.md)           |
| `src/reader/<format>/`                          | one folder per format: `parser.ts` (`parseX`) and its record types                                                 |
| `src/reader/result.ts`                          | `ReaderResult`: `Success` (result, warnings) or `Error` (message, line)                                            |
| `src/reader/common/text/tokenizer.ts`           | `Tokenizer` over a `StringLike`: line, whitespace and field scanning, `Tokens` (start/end pairs)                   |
| `src/reader/common/text/text-byte-tokenizer.ts` | `TextByteTokenizer`: the same functions over UTF-8 bytes, `TextByteTokens`; the text readers' tokenizer            |
| `src/reader/common/text/column/`                | columns over tokens: `TokenColumn`, `FixedColumn` and their byte forms (`TextByteTokenColumn`, …)                  |
| `src/reader/common/text/byte-string-cache.ts`   | `ByteStringCache`: one string per repeated short value of a byte column (cif.md §3)                                |
| `src/reader/common/decode.ts`                   | base64 to bytes (from bytes or a `StringLike`)                                                                     |
| `src/common/binary-cif/`                        | BinaryCIF: `encoding.ts` (file and encoding types), `decoder.ts`, `array-encoder.ts`, `classifier.ts`, `header.ts` |
| `src/common/msgpack/`                           | MessagePack decode (with extensions and skipping) and encode                                                       |
| `src/common/netcdf/`                            | NetCDF reader (NCTRAJ)                                                                                             |
| `src/common/` (rest)                            | `FileHandle` (random-access reads), `SimpleBuffer` and `IOBuffer` (endian-aware), typed-array reads, ASCII helpers |
| `src/writer/cif.ts`, `src/writer/cif/`          | `CifWriter`: categories and fields to text or BinaryCIF (cif.md §6)                                                |
| `src/utils/short-string-pool.ts`                | string interning for short strings                                                                                 |

## Readers: the contract

1. **`parseX(data)` returns a `Task<ReaderResult<T>>`.** Bad input is a `ReaderResult.Error` with a line when the reader
   knows it, not a throw; the task reports progress and can be aborted where the reader chunks (`chunkedSubtask`).
2. **What a reader takes, by kind:**
   - text formats: `StringLike | ReadonlyUint8Array` (CIF, CSV, cube, DX, Dynamo, GRO, LAMMPS, MOL, MOL2, OBJ, PLY,
     PRMTOP, PSF, SDF, TOP, XYZ), tokenized as bytes (§Decisions);
   - text formats that take a `StringLike` only: PDB, XSC, the CryoET Data Portal NDJSON;
   - binary formats: `Uint8Array` (BinaryCIF, DCD, MTZ, NCTRAJ, Simularium, TRR, VTP, XTC, Artiatomi), CCP4 and DSN6
     also a `FileHandle`;
   - RELION STAR: a parsed `CifFile` (STAR is CIF syntax).
3. **A reader only reads its bytes** (`ReadonlyUint8Array`, core `type-helpers`). A parsed file may keep them and read
   values from them on access: a caller that changes the buffer changes the values. The plugin hands readers a
   response's shared bytes without a copy for this reason.
4. **Offsets are byte offsets**, Uint32 in the tokens: text of 4 GiB or more is refused (`TextByteTokenizer.bytesOf`,
   `parseCifText`).
5. **Values are columns.** A reader exposes `Column`s (core `data/db`) or `CifField`s, read per row or as arrays;
   numbers parse from the bytes (`parseIntBytes`, `parseFloatBytes`, core `number-parser`), bit-identical to the string
   parsers. `.` and `?` are CIF's not-present and unknown (`valueKind`), read as `''` and 0.

## Decisions

- **Text is tokenized as bytes.** A downloaded or read file is its bytes; a byte tokenizer neither decodes the whole
  text nor holds a string of it, and indexing bytes is about 2× faster than `charCodeAt` (slower still at a call site
  that has seen a `ChunkedBigString`). The delimiters of every text format here are ASCII, which never occurs inside a
  UTF-8 multi-byte sequence, so scanning bytes is safe; values decode as UTF-8 when read. A string input is encoded once
  up front (`utf8Encode`: one `encodeInto` pass) and takes the same path: one tokenizer per format, not one for strings
  and one for bytes, which would diverge. The cost is the encoded copy while the caller keeps its string.
- **The string `Tokenizer` stays** for code that tokenizes strings (the PDB and NDJSON readers, model's PDB conversion
  in `formats/structure/pdb/`); its byte twin has the same functions, so a reader moves to bytes by changing its import.
- **Decode on access.** Text values are byte ranges until read; BinaryCIF columns decode on first access (cif.md §4). A
  reader of a big file pays for the columns a model reads, not for the file.
- **Fixed-column formats count bytes** (`TextByteFixedColumn`): their columns are ASCII, so bytes are characters.
- **Only CIF is written here** (`CifWriter`); the other formats are read-only in io.

## Tests

`_test/` next to the code (`reader/_test/` holds one file per format). Oracles: the string CIF tokenizer
(`reader/cif/text/_test/string-parser.ts`) for the byte tokenizer's fields, `decodeSteps` for the fused BinaryCIF
decode; fixtures are small and generated in the test, real files only from `data/examples`.
