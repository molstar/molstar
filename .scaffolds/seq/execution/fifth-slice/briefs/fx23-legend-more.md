# Lane FX23: a legend's "+n" opens and closes

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx23`. Model:
Sonnet. Beside you: RX (its fix round: `ext-rondo`, card providers, the ring's hover, bindings), GX (examples' genomes,
`grid-region.ts`, matrix select). Stay out of their files; `viewer/src/legends.ts` and the legend's view are yours.

**The author, 2026-11-23** (a screenshot of the genome ring's legend "Chromosomes (values in view) chr1 … chr12 +18"):
"another issue here: the +18 can't be expanded (and later collapsed)".

**Measure first**: `viewer/src/legends.ts:229-246` caps a legend at `maxItems` (12) and counts the rest into `more`,
drawn as a plain "+n". Also check the count itself: the ring shows GRCh38's 23 chromosomes, yet 12 + 18 = 30 values. Say
which values the legend lists past the 23 (likely the real cytobands' alt contigs or a domain read from a source, not
the view) and whether "values in view" is true. If the legend lists values not drawn, that is a bug here: list only the
drawn ones, with a test failing on the parent.

**Fix**: "+n" is a button (`aria-expanded`, a label "Show n more" / "Show fewer") that opens the legend's full list in
place, wrapping as today, and closes it again. Keyboard (Enter, Space) and pointer alike. The open state is the
session's, per legend (kept across redraws and a pan that keeps the same legend, never written in the spec). A legend
over a very long domain (hundreds) opens to a bounded, scrollable list (a max height), never the whole page. Headless
models first (the legend's model holds `open`, the view renders it; `@molstar/ui` as the other views).

**Tests** (each failing on the parent): the button's label and `aria-expanded`; open lists every entry, close returns to
12 and "+n"; the state survives a redraw; a 500-value legend opens bounded; the count equals the values drawn.

**Where**: a worktree on `seq/fx23`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); one commit (two if
the count is a bug); hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`). No language
change. Anything in `@molstar/core` or `@molstar/io` stops the lane.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers
to completion, `pnpm -s -C packages/seq/core size` (read the rows), `pnpm test:medium packages/seq/`.

**Report**: under 200 words: what was measured (the count too), the fix, commits, checks with numbers, a story id with a
long legend, the `.v6/changelogs/` text.
