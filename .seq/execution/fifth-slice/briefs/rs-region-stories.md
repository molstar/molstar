# Lane RS: region reads' stories (RG's step 5)

**Base**: branch `seq/rg` at `1d638cde5` (RG's work, not yet merged; `git reset --hard seq/rg` first, then branch
`seq/rs`). Model: Sonnet (stories, generators and tests on a built core; the orchestrator's look and a review gate
it). RG's own review runs beside you, read only. No other lane writes.

**Task**: what [region-reads](../../../design/region-reads.md) asks of Storybook and `seqx-examples`
(Consequences: Storybook, `seqx-examples`; Defaults 1–9), on the core RG built (read `data/tables-region.ts`,
`data/summaries-*.ts`, `render/pairs-*.ts`, `render/grid*.ts`, the matrix layer, `describe-region.ts` first):

- **Generators** in `seqx-examples` (small, deterministic on core's `PCG`, reusable by GX's Genomes dataset later): a
  symmetric Hi-C over the toy assemblies with two levels and a one-level twin; a synteny alignment as sided pairs in
  the `seq` format; 12,000 links held and as BEDPE (one picture); a genotype-like matrix of 40 samples with a long
  deletion; values with missing ones. Each a module with its test; no file over a few hundred kB committed (generate at
  test time where you can).
- **Stories** `Features/Region reads/*`, one per Default: pairs zoomed in (items as arcs); pairs zoomed out (bin pairs
  stroked by count, the hover naming the bin pair); a genome ring at the fit (bin-pair chords); a Hi-C map by region
  (tiles, a pending tile drawn coarse); a Hi-C map over the genome at the fit (one read through the whole-assembly
  level); samples × positions (the matrix layer, a long deletion solid); a matrix on a ring; bins with nothing in view
  and the zoom-in state (the hatched band, "zoom in to read this track"). Light and dark, canvas and SVG, a sentence of
  text each saying what is read and drawn. Each with its headless twin on a fake clock (tiles arriving coarse to fine;
  no real-turn waits).
- **A Design page** "Region reads" (`storybook/src/_stories/design/region-reads.mdx`, as the other Design pages: the
  design's intent, decisions and pictures in short, linking the stories).
- **The host-adapter consumer test**: a host adapter (`seq-data-ports`) serving the new shapes (pairs by region, grid
  tiles, a matrix) through the conformance suite and drawn by the viewer, in one test (the two consumers).
- **Matrix state labels** in the header, `order` and `hidden`: only if small and the design says how; else list them.

**Out of scope**: changing core's data or render behaviour (a bug you find there: a test that shows it, then stop and
report it; RG's review will see it too); sizes (the author decides; do not add to the viewer bundles: stories and
generators are not in them).

**Where**: a worktree on `seq/rs`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each green
step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `.seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`): list them. Feature rows for `.seq/design/features.md` and
`design/features.mdx`: list them. No network. Stop any Storybook server you start.

**Checks**: the quick lane for the folders touched (under 2 s each), `pnpm test:unit packages/seq/`, `npx tsc -b
packages/seq/core packages/seq/viewer packages/seq/examples`, each `typecheck:test`, the storybook typecheck,
`pnpm lint`, `pnpm format:check`, the storybook perf tier, `pnpm -s -C packages/seq/core size` (unchanged by you), the
`.seq/_checks` (features-planned with `--plan .seq/execution/fifth-slice/plan.md`). EC ids the stories show
(EC-169, EC-171, EC-208, EC-220 to EC-223) in the twins' titles.

**Report**: under 300 words: commits; story ids; checks with numbers; core bugs found; what is left; hot-file and
features rows; the `.v6/changelogs/` text.
