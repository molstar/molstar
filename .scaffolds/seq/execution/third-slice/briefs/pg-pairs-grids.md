# Lane PG: pairs, grids and the dot plot

**Goal**: build the pairs and grids half of [grids-and-rows](../../../design/grids-and-rows.md), on RW's merged shared
pieces (the y fields of `Hit` and the hover, the targets, the shared reduction `render/pyramid.ts` and the image LRU,
the rows table on 1D row tiles): the `grid` section layout with `strip.y` and `ctx.axisY` (the layouts design's y path),
`of_y`, `follow_y`, `aspect`; the `pairs` annotation and layer (section scope, arcs, chords, dots, split, kinds, hit on the
curve); the `grid` annotation and layer through the shared reduction and images (triangles); the crosshair and the
rectangle selection, rect subtraction in the Location; the query forms (rects over rows × columns, rects subtract).
Plain defaults (`pairs.mirror` false; the tRNA builder sets it). No 2D store this slice (the scope cut).
The examples: tRNA (arcs nested and crossing, the base-pair dot plot with crosshair hover and rectangle selection;
covariation computed on 1EHZ's positions, a grid layer drawing only on its own systems) and the Complex's PAE (one
section per chain: the grid's x follows CDK2's section, its y cyclin A's) and contacts (a `pairs` layer); each three
ways with a headless twin. Its matrix cells (the rectangle, the crosshair), Configuration and Tensions bind.

Not in scope: the rows layer (RW, merged), navigators (NV, beside you), stories (SY).

**Performance**: the design's targets for pairs and grids (a dot plot's frame, the arcs' layout) under `itWithin`, best
of three; draw and hover 0 B warm (add the tier's cases).

**Two consumers**: one grid read by the grid layer and by NV's miniature through the shared reduction; arcs drawn by
the canvas and `renderSvg` (lists equal).
**Read first**: [grids-and-rows](../../../design/grids-and-rows.md) whole; [layouts](../../../design/layouts.md) (the interface, `bandTop`, `rowScale`, the viewport
record, the y path); [ux-foundations](../../../design/ux-foundations.md) section 8; [viewer-models](../../../design/viewer-models.md);
[performance](../../../design/performance.md); [edge-cases](../../../design/edge-cases.md); `packages/seq/AGENTS.md`; the core
and viewer `_spec`s; the code they name.

**Model**: strong (opus). **Base**: `packages-seq` at `51787bb2b`: `git checkout -B seq/pg-pairs-grids 51787bb2b`; check
`git merge-base --is-ancestor 51787bb2b HEAD`. **The language stays 0.4.0** for the whole slice (one bump, unreleased).

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
