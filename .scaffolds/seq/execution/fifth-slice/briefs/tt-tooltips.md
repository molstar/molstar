# Lane TT: tooltips, and the card composed by providers

**Base**: `packages-seq` at the orchestrator's head after HC merges (`git reset --hard packages-seq` first). Branch
`seq/tt`. Model: strong. **The language, alone**: the minor bump to `0.9.0` is this lane's (the version, its history
entry, the catalogs, the builder); no other lane writes in `lang/` while you run.

**Why**: [tooltips](../../../design/tooltips.md), agreed (the author, 2026-11-22: "Yes, as recommended"; the design's
Intent quotes the asks). Read it whole first: Decision, Defaults, Rejected, Edge cases (EC-250 on), Tests,
Consequences, Answered.

**Task**, a green commit at each step:

1. **Card providers**: the hover card composed by `card-provider` entries in the viewer spec, ordered by `order` then
   name as header entries are; the viewer's own sections (place, item, region, pair, chemistry, …) become built-in
   providers writing today's `Parts` (title, rows, words, `whole`) so the live region's words and the card tests do
   not change (a test that every existing card and sentence is byte-equal before and after). A host adds or replaces a
   provider. HC's hover contract holds; a warm hover stays 0 B.
2. **The language**: the `tooltip` node as the design gives it (by data, by a table inline/URI/source/CIF category,
   MVS's `schema`, `field`, `text_format`), the SVG `export.titles` option; the bump to `0.9.0` with its history line;
   the builder's `views.tooltip(...)` and `view.tooltip(...)`; catalogs and config-keys tests. **Not** the matrix's
   `member_labels` (the design's separate change; a later lane).
3. **Resolve and the store**: tooltip tables through the store (inline too), indexed there, their tiles requested with
   the layers they decorate and kept as those are; pending, empty and failed as data-store.md says.
4. **The `seq:tooltips` provider**: answers only on drawn marks (or bins of them), ignores `whole`, one cap shared with
   the place's items (`moreCount`), Markdown inline from `@molstar/core/util/markdown` (parsed once, cached by text;
   cleaned after parsing; `{}` values escaped; only `!seq-…` links, through commands), plain text for the live region,
   headless cards and SVG titles.
5. **Stories**: `Features/Tooltips/*`: by data (a GFF's attributes), by a table of ranges over a protein (inline and
   from a source), over an assembly by name, Markdown with an in-view link, the SVG titles option; light and dark,
   canvas and SVG, headless twins on a fake clock; the builder showcase gains `.tooltip()`.

6. **The card's overflow** (`fixes.md`, "A hover card's long values overflow its box", the author's): values wrap within
   a capped card width, then an ellipsis; a long VCF REF titled by its length, not a hash. A test per case (canvas-free:
   the card's measured box).

**Tests**: the design's Tests section, each EC in a title; the provider order; MVS parity cases (the same annotation
file gives the same texts as MVS's runtime would, by reading); a cold and a warm card show measured and cited.
**Performance**: a warm hover 0 B with tooltips over a million ranges read by region (an allocation case), the
card show's bytes recorded.

**Where**: a worktree on `seq/tt`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each green
step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`): list them. Anything in
`@molstar/core` or `@molstar/io` stops the lane with a report (using `util/markdown` and `util/string-format` as they are
is fine). A file past 800 lines is a question. Sizes soft: trim first, keep the provider machinery lean and tooltip
tables behind a lazy door where rare; report growth with its cause; never edit bounds. Timing and allocation bounds are
the gates.

**Checks**: the quick lane `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, both perf tiers to completion, `pnpm -s -C packages/seq/core size` (read the rows for FAIL),
`pnpm test:medium packages/seq/` (HC's sweep included), the `.scaffolds/seq/_checks` (features-planned with
`--plan .scaffolds/seq/execution/fifth-slice/plan.md`, `defaults-designed.ts`). A fix's test fails on its parent.

**Report**: under 400 words: commits; checks with numbers; story ids; sizes; the version and its history line; designs
to amend (extensibility, extension-openings, security as the design's Consequences list); hot-file needs; the
`.v6/changelogs/` text.
