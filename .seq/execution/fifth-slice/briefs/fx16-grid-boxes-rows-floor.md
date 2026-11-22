# Lane FX16: the grid's selection boxes; the rows' zoom floor and limits

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx16`. Model:
Sonnet. HC (hover emphasis across layers, legends, the allele card) and LA4 (the liftover: maps, bridges) are running:
stay out of hover emphasis (`drawOwnHover`, `tracks-draw.ts`'s hover path), legends, `describe*.ts`, bridges and maps.

**Items** (`fixes.md`, Open; read each whole):

1. **"A selection's boxes on the grid don't meet the cells' edges"**: measure first (the PAE and contacts story, Over
   the complex: a range selected on cyclin A2's pLDDT), then each selection piece clipped to its cell (the blocks' gaps
   cut it), its stroke inside, a row band and a column band's intersection drawn once; canvas and SVG alike. A twin test
   that every piece's box lies inside its cell's box at a block boundary and at the grid's outer edge (fails on the
   parent).
2. **"Zooming a lazy alignment out past a pixel a row leaves a blank plot"**, part (1) only: while a lazy source serves
   no row levels, the rows' zoom-out stops at a device pixel a row (the rail's grip, `[`, the wheel), so the window
   can't reach the blank state (a test that fails on the parent: the rail can't pass the floor on a source without
   levels). Part (2), the coarse row levels, is a design: not here.
3. Same item, **the limits made the host's**: `rows.zoom` in `ViewerConfig` (`minPx`, `maxPx`, `enabled`; defaults
   today's `MinRowScale` 1/64 and `MaxRowScale` 40 from `core/src/render/rail.ts:53`), read by the zoom keys
   (`viewer/src/behaviors/rows.ts`), the rail's grips and its idle rule; `enabled: false` keeps the rail a pan-only bar
   (no grips, zoom keys and wheel-with-modifier inert, scroll kept); the effective floor is the larger of `minPx` and
   the source's (item 2). **No spec override** (that is a language change: listed for later). The core `_spec`'s and the
   viewer config's docs in the same commit.

**Where**: a worktree on `seq/fx16`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); a commit per item;
hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`). No language change. Anything in
`@molstar/core` or `@molstar/io` stops the lane. Sizes soft: report growth; never edit bounds.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers
to completion, `pnpm -s -C packages/seq/core size` (read the rows), `pnpm test:medium packages/seq/`, the `.seq/_checks`
(features-planned with `--plan .seq/execution/fifth-slice/plan.md`). A fix's test fails on its parent.

**Report**: under 250 words: per item what was measured and fixed, commits, checks with numbers, story ids, the
`.v6/changelogs/` text, designs to amend (grids-and-rows.md's zoom rows).
