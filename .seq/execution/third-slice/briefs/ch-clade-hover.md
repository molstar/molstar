# Lane CH: the clade hover and the ruler's doubled label

**Goal**: a hover on a tree node or a summary row lights its clade, as [grids-and-rows](../../../design/grids-and-rows.md)
now says (Clades, "Hovered"; ux-foundations section 8: the element is hovered, here and on linked sections): the
clade's branches in `seq-hover` ink and one band over its rows (folded: the summary row's band), drawn as a header hover
lights a row's band. Found by the author in the Globins example: hovering the tree lights nothing. Cause seen by the
orchestrator: `render/rows.ts` writes a `tree-node` hit with no `position` and no `yPosition` (`writeHit(..., 'tree-node',
cladeHit(...))`), and `drawHover` draws only from those, so a clade hover draws nothing. Check the viewer's hover path too
(the HoverModel and linked sections keep a hit whose element and extent are a clade's).

**Second item, the ruler's doubled label** (the author, from a screenshot of Globins › Rows › Figure, light): in the
wrapped block for columns 61–80 (groups of ten with a gap between), "70" is drawn twice, slightly offset: once raised,
as the end number the "80" also shows, and once as the regular tick label beside "65" and "75". The first group's end
number lands on a tick's label. A position must carry one label: where an end number (a group's or a block's) falls on a
labelled tick, draw one (the end number's place and style win, or the tick's: pick what the rows `_spec` and the axis
design's label rule say, and say which), and no label may overlap another (the ruler's label gap). Test it headless on
that block of the Figure story (one "70" text in the list; no two label boxes overlapping, over every wrapped block of
the Globins figure and at a few widths), failing on its parent. Also check the end numbers sit on the same baseline as
the tick labels (the "61" and "80" look off the "65"/"75" line) against the design; a departure is fixed, a design
silence goes in the report as drift.

Not in scope: anything else in rows; navigators (NV) and pairs and grids (PG) run beside you.

**Tests**: a headless test that hovers a tree node and a summary row of the Globins view (and a small tree fixture) and
checks the drawn list: the clade's branches in the hover ink, one band covering exactly its rows' span, nothing for
another clade; a pinned row of the clade lit at both its places. Shown failing on its parent (name the command). Draw
and hover stay 0 B warm (the allocation tier's rows hover case).

**Model**: fast (sonnet). **Base**: `packages-seq` at `00096372f`: `git checkout -B seq/ch-clade-hover 00096372f`.

**Minimal code**: reuse `hoverBand` and the tree's own branch drawing; no private copies.

**Checks**: `pnpm test:unit packages/seq/`; `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`;
`pnpm -C packages/seq/core typecheck:test`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`.
Green only on exit 0 with a non-zero test count.

**Headers**: a modified file keeps its header, its year range extended to 2026 and `@author David Sehnal
<david.sehnal@gmail.com>` added if missing. No mention of `.seq/` in code; code cites its own `_spec` sections (add the
hover rule to the core `_spec`'s rows section).

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, every `package.json` and `tsconfig*.json`,
anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix, or the fix needs a design change.

**Report**: commits; each check with exit code and test count; the stories to look at (Globins, Family, Folded) and
what to look for, light and dark (the Figure story's ruler: one label a position, no overlap).
