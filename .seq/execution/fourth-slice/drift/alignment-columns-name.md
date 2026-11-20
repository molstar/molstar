# Drift: an alignment's columns system has one fixed name

**What**: lane AC names every alignment's columns system `columns` (`packages/seq/core/src/resolve/model.ts`,
`AlignmentSystems`). A ranges source on an alignment reads its rows by that name. A records source shared by two
alignments (Rfam's Stockholm files hold many families) cannot tell their columns apart. No design names this system
(building-blocks.md, maps.md); the `_spec` line is the lane's decision.

**Options**:
1. Keep `columns`, one alignment a source: a shared source is split by its host. Nothing to build.
2. Name the columns by the alignment's ref (`<alignment>/columns`), so a shared source's records say which. A small
   change in resolve and the ranges source; the language unchanged.
3. A source names its alignment by a param (`:alignment ref`), the system name stays `columns`.

**Recommendation**: 2. Records keep naming what they mean without a param, and one Stockholm file of many families
reads as it is written.

Also owed (not drift, a stopgap): the design's `alignment` map kind (building-blocks.md, the gapped row) is not built;
a chain reaches the columns through a `segments` map computed from the row (the `_spec` says so; TODO.md).
