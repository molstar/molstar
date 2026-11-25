# Lane SC: a range's strand in the hover card

**Goal** (CS's leftover): the hover card and the announcer say a range item's strand when its data has one (a gene on
the minus strand: "Strand −"), for a source's items as for inline annotations (`viewer/src/describe.ts:169-184` shows
title, Extent and Category only; its category is read from inline annotations only). Then the p53 chromosome's index
summary may say "name, strand and span" again (`packages/seq/examples/src/views.ts:140`): the one line in examples you
may touch. Test (fails on the parent): a minus-strand gene from a source by region, its card and its announcement.

**Not in scope**: anything else. CF works in `packages/seq/core` (store, ruler, templates, accessibility, rows-model),
MV in core's renderer and the viewer's composition (`viewer.ts`, `default-spec`, the behaviors' wiring): stay in
`describe.ts`, `behaviors/hover-card.ts`, the announcer and their tests.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/sc-card-strand
<that commit>`.

**Checks**: as [hc-hover-card.md](hc-hover-card.md)'s, without test:perf and size.

**Headers**, **hot files** and **stop when** (a public API change): as every lane's (see
[qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; checks; the test failing on its parent; a `.v6/changelogs/` entry's text.
