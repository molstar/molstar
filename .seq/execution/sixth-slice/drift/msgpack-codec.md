# Drift: where the MessagePack codec lives

**What**: tree-core.md ("Codecs live where their libraries do") puts the MessagePack codec in `@molstar/io`, on io's
`encodeMsgPack`/`decodeMsgPack`, imported by `@molstar/seqx-core` through its allowlist; format.md's Built note said
"in io's `repr` codec". The code has it in `@molstar/core/repr/codecs/msgpack`, on `@molstar/core/data/msgpack`, which
moved down from io; io re-exports both at their old paths (repr spec, "Imports go one way").

**Where**: `packages/core/src/repr/codecs/msgpack/`, `packages/core/src/data/msgpack`; the repr spec's Codecs section
and its folder table say so.

**Options**: (a) the design follows the code: MessagePack is core's, its encoder needing no io dependency; (b) move the
codec back to io.

**Recommendation**: (a). The move took an io dependency away from every reader of JSON's binary twin, and the spec
already says it. tree-core.md's text stays as history, its "In the specs" row pointing here; format.md's Built note was
corrected in this pass.
