# The MessagePack codec of the repr tree

`src/repr/msgpack.ts`: `MsgPackCodec`, a `Codec<Uint8Array>` (core's `repr/codec/codec.ts`; core's repr spec, "Codecs").
It is built on `common/msgpack`'s `encodeMsgPack` and `decodeMsgPack`/`parseMsgPackValue` and their extension hooks;
io's MessagePack is not changed. It lives beside the format code in `src/repr/` (io holds the byte formats; core holds
the tree and the codecs that need no library). Its test runs core's `checkCodec`.

## 1. The bytes

- A node is a map `{ kind, id?, params?, custom?, children?, args? }` (keys in that order), a list an array, a record a
  map without `kind`: the JSON codec's shape, so `json(msgpack(x))` is `json(x)`. A reader tells the encodings apart by
  the first byte (a map, never `{`).
- What MessagePack lacks is an extension (type, payload):

| Type        | What               | Payload                                                                                         |
| ----------- | ------------------ | ----------------------------------------------------------------------------------------------- |
| `0x10`      | NaN, ±Infinity, -0 | one byte: 0 NaN, 1 +Infinity, 2 -Infinity, 3 -0                                                 |
| `0x11-0x18` | typed column       | the elements' bytes, little-endian: int8, uint8, int16, uint16, int32, uint32, float32, float64 |
| `0x19`      | category column    | u32 dictionary size; each entry a u32 byte length and its UTF-8; the int32 codes                |
| `0x1a`      | bigint             | its decimal digits, UTF-8                                                                       |

Why extensions: `encodeMsgPack` writes a non-finite number as nil and -0 as 0, and has no bigint. Why one byte for the
specials: they are rare, the common case is a plain number.

- A column is always typed (`typedColumns: true`): there is no plain form, so the language is not needed to read it.

## 2. Zero copy

- The encoder writes a column's bytes from where they are (a `Uint8Array` view over its buffer): no copy beyond the
  output's. The wire is little-endian; a big-endian host byte-swaps a copy on both sides (`swapBytes`).
- The decoder gives a column as a typed view into its input when the payload's offset is a multiple of the element size
  (always for 1-byte elements), else one copy of the payload. A hit is by the position in the document: the encoder
  cannot align a payload, `common/msgpack`'s hook sees no offset. A tree decoded in place shares its input's memory: the
  caller does not change the input afterwards.
- Decode never writes into its input: the tree makes NaN canonical in a column it owns, so a float column holding a NaN
  is copied (a scan of the view), aligned or not. Why: a foreign encoder may write a NaN payload.

## 3. Reading what is not a tree

`decode` never throws. The bytes that are not one MessagePack value (truncated, malformed, bytes after it) are
`repr.codec.msgpack.syntax`, with the offset reached; only a stack overflow of the parser is `repr.tree.depth_exceeded`
(a read past the end is also a `RangeError`). A key given twice in a map is `repr.codec.msgpack.duplicate_key`, a key
that is not a string `repr.codec.msgpack.key_type` (the decoder keeps the last and stringifies; a second pass over the
document, skipping extension payloads, finds them), each with offset and path. A decoded value that is no tree has the
JSON codec's issues, each with its path: `repr.tree.node.*`, `repr.tree.value.reserved_key`, `repr.tree.value.shape` (a
`bin` is no column: only an extension the decoder made is), `repr.tree.value.tag` (a malformed extension: a column of
the wrong length, a category code out of its dictionary, invalid UTF-8 in a category entry or bigint digits, an unknown
type).

Strings are read as values: a leading U+FEFF is kept (`common/msgpack`'s `str()` reads with `utf8ReadValue`, also for
BinaryCIF), and the extensions' UTF-8 is strict. **Encoding refuses a lone surrogate** (it throws, as an unprintable
node does): UTF-8 cannot hold it and it would not come back.

## 4. Tests

`src/repr/_test/msgpack.test.ts` (quick: conformance, JSON and MessagePack round trip by hash, every column type,
alignment, malformed input, keys, BOM, surrogates) and `msgpack.heavy.test.ts` (decode throughput over the bytes the
document holds: copied columns at least 200 MB/s; a node-heavy document is bound by its many small values, about 45
MB/s, guarded at 20).
