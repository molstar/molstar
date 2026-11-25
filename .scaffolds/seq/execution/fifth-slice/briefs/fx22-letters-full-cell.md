# Lane FX22: a letter centred in its whole cell, not its visible part

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx22`. Model:
Sonnet. Beside you: GX (`seqx-examples/src/genomes/`), RX (`ext-rondo`, `render/pairs-*`, the language), FX19, FX20.
Stay out of their files.

**The author, 2026-11-23** (two screenshots of the p53 sequence row while scrolling, the edge cells' letters squeezed
into their visible slivers): "for the sequence tracks, i think we should always keep the letters inside at the center
of their full box, and not just the visible one. it messes up the spacing and feels really weird when scrolling".

**Cause, found by the orchestrator**: `core/src/render/letters.ts:404` takes a cell's `x0`/`x1` from the axis's
`sx0`/`sx1`, which are cut to the view, so the glyph's centre (`letters-draw.ts:224`, `(x0 + x1) / 2`), its `show`
decision (`letters.ts:432–438`, width against the advance) and the three-letter rule use the visible width. Measure
before changing: confirm in a test that an edge cell's glyph x is the visible part's centre on the parent.

**Fix**: a cell's glyph centred on its whole extent (the axis's unclipped positions), shown by the whole cell's width,
and clipped by the band's box where it overhangs; the fill, hover, selection and conflict boxes as today (cut to the
view). Canvas and SVG alike. Do the same for every row that draws text in cells by the same path (molecule letters,
`text.ts` if it centres in cells): list which you changed and which did not need it.

**Tests** (each failing on the parent): at a pan where the first and last cells are partly out, every glyph's x is its
cell's whole-extent centre, so the spacing between glyph centres is constant across the row (one advance a cell) at
both edges; a cell wider than the view at deep zoom centres its letter off-screen rather than in view; SVG equal to
the canvas draw list. Allocation: the pan's frame bytes unchanged (the perf tier's sequence cases).

**Where**: a worktree on `seq/fx22`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); one commit; hot
files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`). No language change. Anything in
`@molstar/core` or `@molstar/io` stops the lane.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers
to completion, `pnpm -s -C packages/seq/core size` (read the rows), `pnpm test:medium packages/seq/`.

**Report**: under 200 words: what was measured, the fix, rows changed, commits, checks with numbers, a story id to look
at, the `.v6/changelogs/` text.
