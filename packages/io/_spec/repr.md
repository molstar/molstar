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
  output's. The platform is little-endian, as the tree's column hash assumes.
- The decoder gives a column as a typed view into its input when the payload's offset is a multiple of the element size
  (always for 1-byte elements), else one copy of the payload. A hit is by the position in the document: the encoder
  cannot align a payload, `common/msgpack`'s hook sees no offset. A tree decoded in place shares its input's memory: the
  caller does not change the input afterwards, and the tree's own canonicalization of NaN writes there.

## 3. Reading what is not a tree

`decode` never throws. The bytes that are not one MessagePack value (truncated, malformed, bytes after it) are
`repr.codec.msgpack.syntax`; a nesting that overflows the parser's stack is `repr.tree.depth_exceeded`. A decoded value
that is no tree has the JSON codec's issues, each with its path: `repr.tree.node.*`, `repr.tree.value.reserved_key`,
`repr.tree.value.shape` (a `bin` is no column: only an extension the decoder made is), `repr.tree.value.tag` (a
malformed extension: a column of the wrong length, a category code out of its dictionary, an unknown type).

## 4. Tests

`src/repr/_test/msgpack.test.ts` (quick: conformance, JSON and MessagePack round trip by hash, every column type,
alignment, malformed input) and `msgpack.heavy.test.ts` (decode throughput of a large spec of columns, at least 200
MB/s).
