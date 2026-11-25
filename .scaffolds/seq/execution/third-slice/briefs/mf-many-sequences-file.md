# Lane MF: one annotation file, several sequences (a Files story)

**Goal** (the author, 2026-11-19: "could we have annotation for multiple sequences in a same file too?"): it is built
([files.md](../../../design/files.md), "Its systems, by name": a `target` field names each row's system, GFF3's
`seqid`, a TSV's `chrom`; an empty target and a file with no `target` field are on the source's own system; the
records reader, `core/src/data/records.ts`; Complex reads both chains' secondary structure from one file) but no Files
story shows it. Add one view and its story beside the tables (Features › Files): one small TSV fixture (committed,
generated-and-said or trimmed from the p53 and Complex data already here; `LicenseRef-seqx-generated` if generated,
with its provenance entry like the others) whose rows name two sequences by a `target` column (and one row with an
empty target on the source's own, to show the rule), read by one source, two sections each showing only its own rows;
the file shown under the view as the table stories show theirs. The view in `packages/seq/examples/src/views.ts`
(`FilesViews`, the tables' list near line 629) and `files/tables.ts`; the story in
`packages/seq/storybook/src/_stories/features/files/files.stories.ts`.

**Tests** (fail without the view): each section's items are exactly its target's rows; the empty-target row on the
source's own system; the story's headless twin (as the other Files stories have); drawing against the data at the fit.

**Not in scope**: anything else; no core change (a gap found: stop and report). CF works in `packages/seq/core`, MV in
core's renderer and the viewer's composition: stay in examples and storybook.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B
seq/mf-many-sequences <that commit>`.

**Checks**: `pnpm test:unit packages/seq/`, `npx tsc -b packages/seq/examples`, `pnpm -C packages/seq/examples
typecheck:test`, `pnpm -C packages/seq/storybook typecheck`, `pnpm lint`, `pnpm format:check`.

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; the fixture's provenance; a `.v6/changelogs/` entry's text.
