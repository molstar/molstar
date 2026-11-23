# Lane SY: the stories library and the animation port

**Goal**: build [stories](../../../design/stories.md) as agreed: `stateOf`, `plan` and `sample` with the channel kinds in
`seqx-core` (property tests: `sample(plan, 0)` is `from`, `sample(plan, 1)` is `to`, monotone, deterministic,
allocation-free warm); the frame slot over the shown state; the `AnimationPort` with its default (requestAnimationFrame)
and headless (fake clock) adapters — the animation loop is the application's, Seq\* owns none; the view history's jumps
and focus through it (today a jump eases only the first leading section's start and end: every leading section's
viewport, the text anchor, `scale`, `top`, `bandTop`, `rowScale`, a ring's `rotation` with labels hidden and placed once
at its end — the layouts design's story `view` channel, field by field); reduced motion as a cut; a step stored as a
full state or a patch (keyframes and deltas, `story.maxChain` 8), a step that no longer applies shown as such; the p53
story played, scrubbed and exported (an SVG per step), with its headless twin.

Not in scope: rows (RW, beside you: if a channel for `bandTop`/`rowScale` needs RW's rows, build it against the
viewport record's fields and say so), navigators.

**Performance**: a story's frame: `sample` under its target (`itWithin`, best of three), allocation-free warm (the tier).

**Two consumers**: one plan sampled by the viewer (through the port) and by the headless export at the same t, equal.

**Read first**: [stories](../../../design/stories.md) whole; [layouts](../../../design/layouts.md) (the interface, `bandTop`, `rowScale`, the viewport
record, the y path); [ux-foundations](../../../design/ux-foundations.md) section 8; [viewer-models](../../../design/viewer-models.md);
[performance](../../../design/performance.md); [edge-cases](../../../design/edge-cases.md); `packages/seq/AGENTS.md`; the core
and viewer `_spec`s; the code they name.

**Model**: strong (opus). **Base**: `packages-seq` at `43beaa3e3`: `git checkout -B seq/sy-stories 43beaa3e3`; check
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
