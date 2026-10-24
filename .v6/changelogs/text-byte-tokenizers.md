## IO: text formats parsed from bytes

- The line-based text formats read their UTF-8 bytes, not a decoded string: `TextByteTokenizer`
  (`reader/common/text/text-byte-tokenizer.ts`, the `Tokenizer` functions over a `Uint8Array`), `TextByteTokenColumn`,
  `TextByteFixedColumn`, `CifField.ofTextByteTokens`. XYZ, GRO, MOL, SDF, MOL2, PSF, TOP, PRMTOP, LAMMPS data and dump,
  cube, DX, OBJ, PLY, Dynamo tbl and CSV parsers take `StringLike | Uint8Array` (a string is encoded once); the plugin
  hands them a response's bytes (`SO.Data.readBytes`), held once and never decoded whole. Parse and read every column,
  synthetic files (`node scripts/pipeline-bench.mjs parsers`): XYZ 20×100k atoms 370 → 208 ms, GRO 20×100k 392 → 293
  ms, LAMMPS dump 20×100k 392 → 273 ms, LAMMPS data 1M atoms 253 → 170 ms, cube 160³ 226 → 148 ms, OBJ 500k vertices
  379 → 228 ms, PLY 175 → 126 ms; SDF unchanged. Short string values go through a per-column cache shared with text
  CIF (`ByteStringCache`), which turns off for a column whose values do not repeat.
- The byte parsers and readers take `ReadonlyUint8Array` (core `type-helpers`, a `Uint8Array` without its writing
  members, checked by the types only): a parser cannot write into the shared bytes it is given. Core `utf8Read`,
  `utf8ReadValue`, `utf8ReadLong`, `asciiIndexOf`, `asciiStartsWith`, `parseIntBytes` and `parseFloatBytes` accept it.
- [Breaking] CSV follows RFC 4180 (`parseCsv`): fields keep their spaces, doubled quotes read as one, empty fields keep
  the columns aligned, `.` and `?` are values, comments are off by default; see `.v6/plans/breaking-v6-changes.md`.
- `parseDx(data, name, { binary })`: bytes are binary DX unless `binary: false`; a file without an `object 3` record
  fails instead of hanging.
- Fixes: a Dynamo tbl row with too few columns before the first full row misaligned the columns; a PLY format line with
  extra spaces or another version number is read; an ASCII PLY whose last list row has no line break no longer hangs;
  a GRO without decimal coordinates is an error, not zeros.
