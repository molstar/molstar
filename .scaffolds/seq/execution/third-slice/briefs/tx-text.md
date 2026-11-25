# Lane TX: the layout interface, the tracks restated on it, and the text layout

**Goal**: build the first half of [layouts](../../../design/layouts.md), in its "Order of work":
1. **The interface first, as its own commit the orchestrator can merge early**: `render/layout/` with
   `SectionLayoutType`, `Placement`, strips as a struct of arrays, the translate and polar warps with `WarpWriter`;
   `AxisTransform` with groups and a system a strip; **the tracks layout restated on it** (its records unchanged but for
   the `window` header: the restatement's guard); `FrameInput.viewports` with `top` and the viewport record layouts own
   (per-band `bandTop` and `rowScale`, the grid's `y`), follow from the leader's viewport, the key and caches per
   section; `render/window.ts` (the window mark, its hit targets with a `part`) and the `window` pass (the first
   emphasis pass, under the highlights); glyph runs with per-letter y, angle and vertical scale; the layer contract's
   `used`, `scope`, `ctx.span`, `ctx.place`, `ctx.labels`; the y path for `axes: 2` layers (grids-and-rows asks for
   it). Report when this commit is green: RB (blocks and rings) starts on it.
2. **The text layout**: lines, groups of ten, gutters and numbers per row, inline labels, chips, frames, collapse; the
   reflow on resize within its target; ↑↓ stepping lines (and rows, then the next line block, once rows exist); a drag
   selecting in reading order; a drag past a scrolling box's edge scrolling it, clock-timed; the wheel taken by a text
   section only while it overflows; `seq.render.text.lines`; `renderSvg`'s `window` option and the view's optional
   export params; the story `view` channel's fields for a text anchor.
3. **The viewer for both**: `ViewModel` through each section's `SectionNav` (the viewport record's meanings, `scale`,
   `clamp` per layout); `input` reading the placement (`at`, `offsetAt`, the section under a point); the hit target
   `window` with its `part`, labels hit as `item`, gutters as `row`; the `window` behavior provider (its drags, slider
   semantics: `aria-valuetext` on settle, Esc cancels, arrows pan, Shift+arrows resize); the `reveal` group
   (`reveal.cursor` moved from `cursor.reveal`, `reveal.linked`); `selection.copyLetters` copying what is shown
   (reading order, the complement as drawn, shown positions only, rows as FASTA lines; an option for stored letters)
   with its key in `view.keys`; history steps restored by position across a layout change.
4. **The p53 text figure**, three ways, drawing checked against the data, with a headless twin.

The design's Decision (the author's decisions first), its matrix cells, its Configuration table ("whose" copied into
the `_spec`s) and its Tensions bind; Tensions stay tensions.

Not in scope: rings and blocks (RB, after your interface commit merges), the rows layer (RW), the ESPript figure (after
RW), navigators (NV).

**Read first**: [layouts](../../../design/layouts.md) whole; [axis](../../../design/axis.md) (what reads the axis; a
break zero width, a hairline, on every layout); [layout-and-draw-lists](../../../design/layout-and-draw-lists.md);
[viewer-models](../../../design/viewer-models.md); [ux-foundations](../../../design/ux-foundations.md) section 8;
[grids-and-rows](../../../design/grids-and-rows.md) (what it asks of layouts: the y path, `bandTop`, `rowScale`);
[stories](../../../design/stories.md) (the `view` channel); [performance](../../../design/performance.md);
`packages/seq/AGENTS.md`; the core and viewer `_spec`s; `packages/seq/core/src/render/` (tracks.ts is the hot file).

**Model**: strong (opus).

**Base**: `packages-seq` at `61d2887c3`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/tx-text 61d2887c3`, and check `git merge-base --is-ancestor 61d2887c3 HEAD`.

**Touches**: `packages/seq/{core,viewer,examples}/src/**`, `packages/seq/storybook/src/_stories/**` (not the Design
pages), the `_spec`s, `packages/seq/TODO.md`, `.scaffolds/seq/design/performance.md`'s rows only through your report.

**Two consumers**: one text section drawn by the canvas and by `renderSvg` (Canvas2D and SVG lists equal); the tracks
records unchanged by the restatement (both before and after).

**Tests**: the design's Consequences, Tests, all of it for the interface, tracks and text (property tests over random
sections a layout; an unchanged section not laid out again; drawing against the data at the fit view, across a tile
edge and across a line break; reused buffers against fresh ones; the p53 text figure's record; the tracks records as
the guard; the nav's arithmetic in Node; the matrix cells as binding tests; back across a layout patch).

**Performance representative**: text reflow under 8 ms a frame for the p53 figure and a 35,000-residue sequence, flat in
length (35,000 against 3,500 at most 1.2×); hover under 2 ms; each under `itWithin` (best of three); draw and hover
0 bytes warm through the warp (the opt-in tier: run it and report).

**Whose each choice is**: the design's table; nothing configurable beyond it; plain defaults.

**Minimal code**: one interface, the tracks restated rather than duplicated; reuse the axis transform, tiles, hits and
the hover index; no private copies; sizes reported (guides, not bounds).

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
design, or a target cannot be met without changing a design (report the measurement and the options), or a permission
is refused (report the action and why; do not work around it).

**Report**: the commits (name the interface commit); each check with its exit code and test count; the sizes and the
performance numbers; the hot-file changes you need; the changelog entry's text; config items with their defaults; the
matrix cells filled; drift (what, where, options, your recommendation); followups. You cannot see a browser: name the
stories to look at, and what to look for in each (light and dark).
