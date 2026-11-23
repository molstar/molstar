# Lane RB: blocks and rings

**Goal**: build the second half of [layouts](../../../design/layouts.md) on the layout interface TX merged first:
- **Rings** ("Rings"): sectors on the polar warp, labels with leaders, a feature across the origin by its shortest arc,
  a break zero width with a hairline (axis), the ring's window and its hit targets through `render/window.ts` (the
  `window` target with a `part`, the shared window behavior: slider semantics, `aria-valuetext` on settle, Esc
  cancels, arrows pan, Shift+arrows resize), a linear detail following the ring's window and a ring following a detail
  (EC-108, EC-114), the cyclic wrap for leaders, `rotation` as a story channel tweened with labels hidden and placed
  once at its end, `seq.render.ring.image`; labels placed only when the ring is laid out.
- **Blocks** ("Blocks"): `block`, `arrange`, `copies` (its plain default with the design's reason), a path a block; in
  all three layouts (tracks, text, rings); `seq.render.follow.blocks`.
- The viewer for both: `SectionNav` per layout (the viewport record's meanings, `clamp`, the cyclic wrap), input
  through the placement, history steps restored by position across a layout change.
- **The Plasmid ring** with its linear detail and the rotated linear view, and the Complex's chains as blocks with the
  ligands as chips (the datasets in [examples](../../../design/examples.md); their fixtures may be minimal stand-ins
  until the dataset lanes, named as such), three ways, drawing checked against the data (across the origin too), with
  headless twins.

The design's Decision (the author's decisions first), its matrix cells (the ring's window), its Configuration table
("whose" copied into the `_spec`s) and its Tensions bind; Tensions stay tensions.

Not in scope: the text layout (TX, beside you: do not edit its files; where you both touch the viewer's nav or input,
keep your edits to rings and blocks, and the orchestrator merges TX before you, so rebase onto `packages-seq` when told),
navigators (NV), the rows layer.

**The language version stays 0.4.0** for the whole slice (one bump, unreleased): add `block`, `arrange`, `copies` and
the rings' params under 0.4.0. **The allocation tier**: "p53 on its gene, whole protein" is over its 64 KB bound at the
base already (a stage 1 regression being traced): do not move that bound; report against the base's numbers.

**Read first**: [layouts](../../../design/layouts.md) whole; [axis](../../../design/axis.md);
[navigators](../../../design/navigators.md) (the window's shared rules: grips, floor, idle, hidden);
[viewer-models](../../../design/viewer-models.md); [ux-foundations](../../../design/ux-foundations.md) section 8;
[stories](../../../design/stories.md) (the `view` channel, `rotation`); [examples](../../../design/examples.md)
(Plasmid, Complex); [edge-cases](../../../design/edge-cases.md) (EC-108, EC-114); `packages/seq/AGENTS.md`; the core
and viewer `_spec`s; `packages/seq/core/src/render/layout/` (the interface).

**Model**: strong (opus).

**Base**: `packages-seq` at `9cfc6801f`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/rb-blocks-rings 9cfc6801f`, and check `git merge-base --is-ancestor 9cfc6801f HEAD`.

**Touches**: `packages/seq/{core,viewer,examples}/src/**` (rings and blocks), `packages/seq/storybook/src/_stories/**`
(not the Design pages), the `_spec`s, `packages/seq/TODO.md`.

**Two consumers**: a ring and its linear detail on one store, each following the other; the ring drawn by the canvas and
by `renderSvg` (lists equal).

**Tests**: the design's Consequences, Tests, for rings and blocks (property tests over random sections: every position
placed once, `at` the inverse of `place`, `offsetAt` monotone in angle; a detail following a ring's window and back; an
unchanged section not laid out again; drawing against the data at the fit view, across a tile edge and across the
origin; reused buffers against fresh ones; the Plasmid ring's record; the nav's arithmetic in Node; the matrix cells
as binding tests; back across a layout patch).

**Performance representative**: a ring's frame (the Plasmid, four rings, about 30 labels: a pan of its linear detail,
the ring drawn from kept lays) under 8 ms; hover under 2 ms; each under `itWithin` (best of three); draw and hover
0 bytes warm through the polar warp (the opt-in tier: run it and report).

**Whose each choice is**: the design's table; nothing configurable beyond it; plain defaults.

**Minimal code**: build on the interface and the shared window; reuse the axis transform, tiles, hits, the hover index,
the label placement; no private copies; sizes reported (guides, not bounds).

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `pnpm test:heavy packages/seq/` (before it, `pgrep -fl "vitest.*heavy"`: wait if another is running)
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`, and `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`, `pnpm format:check`
- `pnpm -C packages/seq/core size` (report the numbers)

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.scaffolds/seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report);
- anything in `.scaffolds/seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
design, or the interface must change in a way TX's text layout would feel (report it: the orchestrator coordinates), or
a target cannot be met without changing a design (report the measurement and the options), or a permission is refused
(report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the sizes and the performance numbers; the
hot-file changes you need; the changelog entry's text; config items with their defaults; the matrix cells filled;
drift (what, where, options, your recommendation); followups. You cannot see a browser: name the stories to look at,
and what to look for in each (light and dark).
