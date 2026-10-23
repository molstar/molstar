# 13: LAMMPS and large text files (notes, not audited)

Notes for a later pass on the LAMMPS readers (`packages/io/src/reader/lammps`), written alongside the CIF parser pass
(byte tokenizer for text CIF, lazy and fused BinaryCIF decoding; see [cif-parsing](../designs/cif-parsing.md)). The
LAMMPS code is unchanged. Apple M4, Node 25, 2026-10-23.

## Measurements

Scanning a text body with numbers, 76 M chars (a LAMMPS-like `id mol type x y z` dump, 2M lines; tokenize and parse the
last three columns per line as floats):

| Input                                                                           |   Time |
| ------------------------------------------------------------------------------- | -----: |
| primitive string (`charCodeAt`)                                                 | 265 ms |
| `ChunkedBigString`                                                              | 511 ms |
| primitive string after `ChunkedBigString` was seen at the same call site        | 363 ms |
| `Uint8Array` scanned directly, inline number parser (not the exact `parseFloat`) |  84 ms |
| `TextDecoder` on ASCII                                                          |   7 ms |

The CIF pass measured the same effect on the real tokenizer: flat-192-cif (124 MB) parses in 285–327 ms over bytes
against 594–726 ms over a string, with the exact number parsers on bytes (`parseFloatBytes`).

## A missing guard: 2³² chars

`ChunkedBigString` maps an index to its chunk by `index >>> STRING_CHUNK_SHIFT` (`core/util/string-like.ts`), and
`TokenBuilder` stores offsets in a `Uint32Array` (`io/reader/common/text/tokenizer.ts`). Past 2³² chars both wrap around
without an error: a file of 4 GiB or more reads garbage. The CIF byte tokenizer refuses such input; the LAMMPS readers
need the same guard (or Float64 offsets).

## `traj/parser.ts` `handleAtoms`

- Allocates 6–9 token arrays per frame (`id`, `mol`, `type`, `x`, `y`, `z`, more with velocities), each a start and an
  end per atom: 48–72 B per atom-frame against about 40 B of text. All frames are tokenized before any is converted.
- `mol` and `type` tokens are needed only for frame 0.
- The column `switch` (`case columnIndexMap['id']:` ...) re-reads `columnIndexMap[...]` properties per token.

Example: 100k atoms × 500 frames peaks at about 2.5 GB of text + 2.4–3.6 GB of tokens + 0.6–1.2 GB of coordinates.

## Options, in order

1. **Parse straight into typed arrays** (coordinates per frame as Float32, ids and types once), no tokens; plus the
   4 GiB guard. Removes the token peak; the text stays.
2. **Lazy frames over bytes**: one pass builds a frame index (byte offset per frame), `getFrameAtIndex` parses a frame
   as a task when asked. `SimulariumParticleTrajectory` is the precedent. Memory: the bytes and the frames in use.
3. **File-backed frames** for files larger than memory: a Data object holding a `File`, frames read and parsed
   asynchronously by offset.
4. **The generic byte tokenizer** that this pass starts for CIF (`reader/cif/text/field.ts`, `parseFloatBytes`), shared
   by the line-based text formats.

## Later: a byte-backed `StringLike`

An ASCII-only `StringLike` over the file's `Uint8Array` (`charCodeAt(i)` = `bytes[i]`, `substring` decoding a range),
opt-in per format provider, so opened files stay `Data.String` while their text is never decoded to a JS string; the CIF
byte tokenizer would detect it and read the bytes directly, skipping `utf8Encode`'s copy. Risks: code that assumes a
small file is a primitive string, and a second `StringLike` kind at the other parsers' `charCodeAt` call sites (a
polymorphic site is ~1.4–2× slower, as `ChunkedBigString` showed above).

