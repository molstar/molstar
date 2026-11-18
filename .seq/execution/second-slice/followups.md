# Followups

- (M) `maxBytes` on `file://` reads in Node (`fetchData`'s file path ignores the four options).
- (M) The asset conformance suite sees a manager through `fetch`; a manager that downloads by other means is outside it.
- (L) The real p53 fixtures round-tripped JSON ↔ MessagePack in a seq package test (io cannot import them).
- (L) `repr.codec.msgpack.syntax` cataloged in core's `CodecIssues` beside the JSON and S-expression codes.
- (L) Duplicate map keys in MessagePack keep the last silently (JSON reports them): an option in `common/msgpack` later.
- (L) Payload alignment is by position: pad in `common/msgpack` if zero-copy decoding ever matters more.
