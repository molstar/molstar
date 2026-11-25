# Lane DP: seq data ports built

**Task**: build [seq-data-ports](../../../design/seq-data-ports.md) (Decision: the port is the source protocol, Seq io, What a source serves for a structure, A host's updates, Annotations given back, Layering; Consequences): `OpenedSource.changes?` and `RegionRequest.field?` (both optional, every runtime unchanged), the store's `reattach` (node content replaced by key, identity kept), `withAdapter` in `context/` (one adapter a format, conflicts listed; resolver and overlay bind by format, no order), `io/structure` as a new `seqx-core` entry (`StructureDescription`, `structureNodes`, keys independent of row order), the viewer's store model reopening a source in place a frame after `changes`. The code keeps `FormatRuntime`; the designs call it the default adapter. No package of its own, no Mol\* adapter, no mmCIF reading (MC), no fake host (HD). The two optional protocol members are public API: list the surface and stop for the author's acceptance before merging.

**Output**: the changes above; one conformance suite every adapter passes (the protocol's plus: equal nodes from a description and a hand-built one, keys unchanged when a chain is removed, `reattach` keeps identity, `changes` coalesced on a fake clock, a malformed host node an issue); properties from the design's Consequences; `_spec` Data (Protocol, Store), Context (Registry entries), a Seq io section, viewer StoreModel. The round trip through `structureFacts` waits for ST1.

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/dp-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
