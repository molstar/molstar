# Lane TE: an item from an earlier tile, lost past its edge

**Goal**: lane TV found that a range starting in an earlier tile is not drawn when the view lies wholly past that tile's
edge: on the synthetic 3M source, a tracks view over [1040, 1080] lacks `r1000`, which spans it (rings and text too).
The synthetic source answers only items that start in a region. An item in view with no mark is a lost item (the axis
design; "what is drawn is checked against the data"). Decide where the fix belongs (the source contract: a region's
answer includes items overlapping it, or the reader: it reads back far enough for the longest item, as the store's
want sets allow), following the data-store and protocol designs (`.scaffolds/seq/design/data-store.md`, `protocol.md`,
`layout-and-draw-lists.md`, the `render/ranges.ts` dedupe by key and start from the second slice's WRAP53 fix), and
fix it there with the least code. Check the real sources (the p53 region, inline tables, any tiled reader) for the
same miss.

**Test first**: a quick-lane test drawing against the data that fails on the base: a view wholly past a tile edge
inside a long item (tracks, text, and a ring view over the same source), every item overlapping the view has its marks,
drawn once. Then the fix.

**Model**: fast (sonnet). **Base**: `packages-seq` at `6ee30a042`: `git checkout -B seq/te-tile-edge-items 6ee30a042`; check
`git merge-base --is-ancestor 6ee30a042 HEAD`. Lane RB runs beside you in `render/` (rings, blocks): keep your edits to the
data reading and the source, small and local; if the fix must touch layout code, say so in your report.

**Checks**: `pnpm test:unit packages/seq/`; `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`;
`pnpm -C packages/seq/core typecheck:test` (and viewer, examples); `pnpm lint`; `pnpm format:check`;
`pnpm -C packages/seq/core test:perf` (no scenario worse). No heavy run while `pgrep -fl "vitest.*heavy"` shows one.

**Headers**: a modified file keeps its header, year 2026, `@author David Sehnal <david.sehnal@gmail.com>` added if
missing; a new file gets the Mol\* header. No mention of `packages/sequences` or `.scaffolds/seq/` in code. **Hot files**: never
stage `package.json`, `tsconfig*`, the lockfile, `.v6/changelogs/` (write the entry's text in the report), `.scaffolds/seq/`.
Update the `_spec` section the fix changes.

**Report** (under 300 words): where the fix belongs and why, the commits, the test failing on the base (the command),
each check with its exit code and test count, the tier's numbers for any scenario that moved.
