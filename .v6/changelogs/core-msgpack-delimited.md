## Core: MessagePack and delimited records

- MessagePack moved down from `@molstar/io` to `@molstar/core`: `@molstar/core/data/msgpack` (`encodeMsgPack`,
  `decodeMsgPack`, their extension hooks) and `@molstar/core/repr/codecs/msgpack` (`MsgPackCodec`, the tree's codec,
  and `msgPackColumnDecoder`, its column extensions as a decoder hook for a plain value). `@molstar/io` re-exports both
  at their old paths (`common/msgpack/decode`, `common/msgpack/encode`, `repr/msgpack`).
- `@molstar/core/data/text/delimited` (`readDelimited`): RFC 4180 delimited records with plain string cells (a
  one-character delimiter, quote and comment, a BOM skipped, LF/CRLF/CR, errors with their line). io's CSV reader is
  unchanged.
