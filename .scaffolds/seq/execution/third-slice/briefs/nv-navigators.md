# Lane NV: navigators

**Goal**: build [navigators](../../../design/navigators.md) as agreed, on the merged layouts (the shared `window`
behavior, `render/window.ts`, `SectionNav`) and RW's rows (the shared reduction and image LRU, the rows' summaries): its
Decision (the author's decisions first: the ruler rail is the default navigator of a linear section, with its grips,
floor, idle, grown and hidden cases, the grab zone's split on an alignment's ruler; a ruler drag pans), all four forms
of EC-110 (the strip; the miniature; the ideogram with the `bands` layer and its `giemsa` scheme; the vertical overview
of an alignment's rows, `axes: rows`) and the ruler rail and slim rail; its matrix cells (the navigator row), its
Configuration and its Consequences (`seqx-core`, `seqx-viewer`, the examples: p53's gene region with chr17's ideogram,
p53's tracks with a strip, the protein's strip showing its follower's window, the Globins 10,000-row view with its
vertical overview and a miniature, the plasmid's linear detail with a strip wrapping the origin; the stories; the
tests). Export choices: `navigators` an optional spec param of a view that `renderSvg`'s option overrides.

Not in scope: rows (merged), pairs and grids (PG, beside you: share the reduction, do not copy it), stories (SY).

**Performance**: the design's: a window move redraws in 0.5 ms and allocates nothing warm; a strip of p53's tracks in
4 ms, chr17's ideogram in 1 ms, the 10,000-row overview from its summary grid in 8 ms; each under `itWithin`.

**Two consumers**: a navigator and the view it moves (drag the window, the view moves; move the view, the window
follows), headless on a fake clock; the navigator drawn by the canvas and by `renderSvg` with `navigators`.
**Read first**: [navigators](../../../design/navigators.md) whole; [layouts](../../../design/layouts.md) (the interface, `bandTop`, `rowScale`, the viewport
record, the y path); [ux-foundations](../../../design/ux-foundations.md) section 8; [viewer-models](../../../design/viewer-models.md);
[performance](../../../design/performance.md); [edge-cases](../../../design/edge-cases.md); `packages/seq/AGENTS.md`; the core
and viewer `_spec`s; the code they name.

**Model**: strong (opus). **Base**: `packages-seq` at `51787bb2b`: `git checkout -B seq/nv-navigators 51787bb2b`; check
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
