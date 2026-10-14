# Vendored: Mol*'s MessagePack

This copy is the old stack's `msgpack` folder, unchanged; the old copy was
deleted at the switch-over (SW1).

`encode.ts`, `decode.ts` and `utf8.ts` are Mol*'s MessagePack codec, vendored so that SeqViewSpec documents and data
can be read and written as MessagePack (seqspec:format#values) without a dependency.
The folder is meant to go back into Mol* with SVS: it imports nothing outside itself.

- **Origin**: Mol* `src/mol-io/common/msgpack/encode.ts`, `decode.ts` and `src/mol-io/common/utf8.ts`, at upstream commit
  `1ebe24097` (Mol* 5.11.0).
- **License**: MIT, Copyright (c) 2017 - now, Mol* contributors; the text in `LICENSE`. The files keep their headers
  (adapted from mmtf-javascript, Alexander Rose and David Sehnal; the decoder loosely based on msgpack-js, MIT, Tim
  Caswell).
- **Ours**: `index.ts`, the typed-array extension, which may go to Mol* as well.
- **Formatting**: Mol*'s (4 spaces, single quotes), so a diff against upstream shows only the changes below; Biome
  neither formats nor lints these three files (an override in `biome.json`).

## Changes

Each change carries a `STORIES-UNTOLD(Mn)` comment.

- **M1, imports and the UTF-8 helper**: `'../utf8'` became `'./utf8.ts'` (Node's type stripping); `utf8.ts` leaves
  out `utf8ReadLong` and so its import of `./string-like`, which the codec does not use.
- **M2, extension types** (Mol*'s codec had none). The encoder takes an optional `MsgPackExtEncoder`, asked for every
  object and array before it is encoded; a `{ type, data }` it returns is written as fixext 1, 2, 4, 8 or 16 for those
  payload lengths, else ext 8, 16 or 32. It is asked twice per value (size, then write), so it should be cheap. The
  decoder reads all six forms (`decodeMsgPack(bytes, ext)`, kept on `MsgPackState.ext`) and passes the type and a
  view of the payload to an optional `MsgPackExtDecoder`; without one, or when it returns undefined, the value is
  `{ type, data }` with a copy of the payload. `skipMsgPackValue` skips them.
- **M3, JSON's rules for what JSON cannot hold** (a document reads the same from either form). Mol*
  threw for a number not finite and for an integer outside the 32-bit ranges, and wrote an undefined object value as
  nil. Now a number not finite is nil (JSON.stringify's `null`), an integer outside int32 and uint32 is a float 64
  (exact, as a JSON number is a double), and an object's keys whose value is undefined are left out.
- **M4, characters outside the Basic Multilingual Plane** (`utf8.ts`): `utf8Write` and `utf8ByteCount` read UTF-16
  code units, so a surrogate pair (an emoji in a description) became two 3-byte sequences, which the decoder reads as
  replacement characters. They read code points (`codePointAt`, the pair's second unit skipped), as the existing
  4-byte branch intended. A lone surrogate is still written as 3 bytes and read back as U+FFFD. Worth fixing upstream.

- **M5, untrusted lengths** (2026-10-14, the review gate's GI-03 and GI-10; worth fixing upstream): `parse` throws
  "truncated" at the end of the buffer (the byte past it read as a positive fixint, so a cut file decoded as data),
  and `array`, `map`, `bin`, `str` and `ext` throw when their length header claims more items than the bytes left can
  hold (an 8-byte document claiming 2^32 elements killed the process or froze it for minutes). An array grows as its
  items are read rather than being made at its claimed length (nested headers that each claimed every byte left
  multiplied by the depth: 1 MB ran the heap out of memory; the lane A review's LA-01). A map's `__proto__`
  key is an own property, as `JSON.parse` makes it, not the object's prototype.

Typed arrays: Mol*'s encoder writes a `Uint8Array` as `bin` and any other typed array as a map of its indices (it is
an object, not an array). With `index.ts`'s `encodeTypedArrays` as the extension encoder, every typed array is the
typed-array extension instead (a `Uint8Array` included, since `bin` is not a JSON type); the decoder still reads `bin`
as a `Uint8Array`.
