# Lane RW: rows, alignments, trees and the shared pieces

**Goal**: build the rows half of [grids-and-rows](../../../design/grids-and-rows.md), first in its "Order of work"
(it lands the shared changes PG and NV build on): its Decision (the author's decisions first: clade folds are session
state with `clades.keep`, a node click toggles a fold, a reveal opens a folded clade; the scope cut: `cells` on 1D row
tiles, no 2D store this slice), its matrix cells, its Configuration, and its Consequences' **Rows** and the rows parts of
`seqx-core`, `seqx-viewer`, the examples (Globins: the family with logo and tree, folded clades, the 10,000-row view;
the ESPript figure's rows in text, with the text layout merged), the stories and the tests. In short: `Hit`'s and the
hover's `ySystem`/`yPosition` and the targets (`cell` with them, `row`, `tree-node`, `summary-row`); the shared
reduction `render/pyramid.ts` and the image LRU on core's `LRUCache` (reused by PG's grid and NV's navigators: make
them general); the column measures; the aligned FASTA, Stockholm and Newick readers with fuzz tests; the `alignment`
and `tree` kinds; the rows layer (order, folds, pins, header, reduction, lazy rows for 10,000 × 300 with the store's
want sets and bounded layout); the logo and the profile; the clade and rows commands with keys in `view.keys` (row
zoom `]`/`[`); the cursor's `part` and `y`; the selection's rows part; text with rows (↑↓ through a line block's rows,
TV's seam); history's record per rows band.

Not in scope: pairs, grids and the grid section (PG, after you), navigators (NV), stories' library (SY).

**Performance**: the design's "Performance and size" targets for rows (the 10,000 × 300 alignment: a pan frame and a
row-scroll frame under their targets with `itWithin`, best of three; layout allocation bounded; draw and hover 0 B
warm).

**Two consumers**: an alignment drawn by the canvas and `renderSvg` from one spec (lists equal, folds passed to the
export); the rows' shared reduction read by the rows layer and a stand-in second reader (the navigators' will be NV's).

**Read first**: [grids-and-rows](../../../design/grids-and-rows.md) whole; [layouts](../../../design/layouts.md) (the interface, `bandTop`, `rowScale`, the viewport
record, the y path); [ux-foundations](../../../design/ux-foundations.md) section 8; [viewer-models](../../../design/viewer-models.md);
[performance](../../../design/performance.md); [edge-cases](../../../design/edge-cases.md); `packages/seq/AGENTS.md`; the core
and viewer `_spec`s; the code they name.

**Model**: strong (opus). **Base**: `packages-seq` at `43beaa3e3`: `git checkout -B seq/rw-rows 43beaa3e3`; check
`git merge-base --is-ancestor 43beaa3e3 HEAD`. **The language stays 0.4.0** for the whole slice (one bump, unreleased).

**Whose each choice is**: the design's Configuration table, copied into the `_spec`s; nothing configurable beyond it;
plain defaults, never inferred from what a section draws. **Minimal code**: reuse what exists (the layout interface,
the shared window behavior, `SectionNav`, tiles, hits, `@molstar/core`'s `LRUCache`, `PCG`); no private copies; sizes
reported (guides, not bounds).

**Tests**: the design's Consequences, Tests list, all of it; drawing against the data at the fit view and across a tile
edge; headless twins on a fake clock for every story; the matrix cells as binding tests; each fix with a test shown
failing on its parent (name the command).

**Checks**: each its own command, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/`;
`pnpm test:heavy packages/seq/` (before it, `pgrep -fl "vitest.*heavy"`: wait if another runs); `npx tsc -b
packages/seq/core packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples>
typecheck:test`; `pnpm -C packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C
packages/seq/core test:perf` (every scenario passes; add the design's allocation cases); `pnpm -C packages/seq/core size`.

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file
for more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.scaffolds/seq/` in code,
comments, tests or `_spec`s; code cites its own `_spec` sections.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`; `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report);
every `package.json` and `tsconfig*.json` (ask in the report: a new entry, `_test` export nulls); anything in `.scaffolds/seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
design, or a target cannot be met without changing a design (report the measurement and the options), or a permission
is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the sizes and the performance numbers; the
hot-file changes you need; the changelog entry's text; config items with their defaults; the matrix cells filled;
drift (what, where, options, your recommendation); followups. You cannot see a browser: name the stories to look at,
and what to look for in each (light and dark).
