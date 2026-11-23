# Drift: unknown kinds of a newer version

**Already recorded**: `packages/seq/TODO.md`, under `repr`, lists `fallback` handling, a version-migration hook and a
version-aware validate ("an unknown kind in a newer minor version is a warning"). This note only ties format.md's
decision to that entry; nothing new is owed.

**What**: format.md ("Unknown kinds") says a kind of a newer minor version, like one in an extension namespace the
reader lacks, is kept, reported and skipped with the kind's declared fallback, never a failure of the whole spec. The
code keeps and skips only a namespaced kind (`repr.schema.kind.extension`, a warning); an unknown kind without a
namespace is an error, and a kind's `fallback` is documented in the description but not acted on (repr spec,
Validation and Schema). In Seq\*, a newer *layer type* is read as its fallback (core spec, Language, Extensions), but a
newer *kind* fails `checkSeqSpec`.

**Where**: `packages/core/src/repr/schema/` (validation of kinds); `packages/seq/core/src/lang/check.ts` (the gate
order).

**Options**: (a) build the TODO entry: an unknown bare kind in a spec of a newer minor of the same major is a warning,
kept and skipped, its `fallback` used where the language declares one; (b) keep the code until 1.0: pre-release `0.x`
readers refuse what they do not know; (c) narrow the design to layer types and extension kinds.

**Recommendation**: (b) now, (a) before 1.0. While the language is `0.x` every minor may break, so a reader cannot
promise a newer minor reads; once it is 1.x the design's rule is what lets an older viewer open a newer figure. The TODO
entry stands as it is.
