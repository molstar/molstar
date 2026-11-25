## IO: MessagePack

- `@molstar/io`: a MessagePack codec for the repr tree (`MsgPackCodec`, `@molstar/io/repr/msgpack`). Typed columns are
  written as extension types 0x11-0x18 and decoded in place where aligned; it also covers category columns, NaN, the
  infinities, -0 and bigints, and round-trips with the JSON codec by hash. Duplicate and non-string keys are reported,
  and a lone surrogate is refused at encode.
- MessagePack strings and keys now keep a leading U+FEFF (BinaryCIF included).
