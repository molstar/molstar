# Lane AX: the axis (direction, collapsed introns, letters through a path)

**Goal**: build [axis](../../../design/axis.md) as agreed, all of it: its Decision (the author's decisions first: the
direction is the writer's, plain defaults, `pieces` `whole` for every section, a break zero width on every layout),
its matrix cells, its Configuration table (the "whose" column copied into the `_spec`s), and its Consequences, which
are this lane's checklist. Its Tensions stay tensions: do not resolve one in code; a drift note if the code forces it.

In short:
- the language: `section.direction` (`genomic`, `reversed`, `follow`, `{ along }`), `section.pieces`,
  `section.collapsed`, `section.pad`; `letters`' `complement` and `code`; the checks;
- resolve: fine-placed points (each nucleotide of a split codon), the coding runs of a codon map (`projectPoint` is
  already exported by C0);
- view: the axis's direction and pieces, its breaks with their hidden ranges, the issues;
- render: `AxisTransform` both ways (extended, not duplicated), breaks and their ticks, letters through one-to-one and
  codon paths, ranges mirrored, the direction mark, residue and break hits, the SVG's data attributes, `renderSvg`'s
  viewport in positions;
- the viewer: offsets and positions through the axis, follow and lead, history steps restored by position (clamped to
  what is drawn, never silently dropped), the cursor and selection across breaks, focus into a hidden stretch and its
  announcement (the cursor names a hidden stretch it crosses), the `break` hit target, `section.reverse` and
  `section.collapse` (patches) with their keys in `view.keys`;
- the examples: each p53 view spells its direction and pieces in its builder (the application drives the
  composition): the gene following its protein, the protein and its gene on one axis (`{ along: p53 }`,
  `pieces: coding`), the genomic follower for contrast, the reversed gene with its complement; each three ways, with
  headless twins.

Not in scope: themes and colors beyond the break ink token's use, layouts other than tracks (a ring's break is RB's;
the rule "zero width, a hairline" is shared), navigators, anything in `packages/core`.

**Read first**: [axis](../../../design/axis.md) whole; [maps](../../../design/maps.md),
[viewer-models](../../../design/viewer-models.md), [layout-and-draw-lists](../../../design/layout-and-draw-lists.md),
[ux-foundations](../../../design/ux-foundations.md) section 8, [edge-cases](../../../design/edge-cases.md) (EC-31,
EC-32, EC-33, EC-117, EC-144, EC-233); `packages/seq/AGENTS.md`; the core and viewer `_spec`s; the code they name.

**Model**: strong (opus).

**Base**: `packages-seq` at `162c745be`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/ax-axis 162c745be`, and check `git merge-base --is-ancestor 162c745be HEAD`.

**Touches**: `packages/seq/{core,viewer,examples}/src/**`, `packages/seq/storybook/src/_stories/**` (the stories; not
the Design pages), the `_spec`s, `packages/seq/TODO.md`.

**Two consumers**: the protein and its gene on one axis drawn by the canvas and by `renderSvg` from one spec, their hits
and the export's `data-*` agreeing; a follower on a reversed collapsed axis moved by its leader and leading it back.

**Tests**: the design's Consequences, Tests list, all of it (property tests, the direction param, hit-testing through
each reader, the view history, drawing against the data, headless twins and records).

**Performance representative**: a pan frame over a reversed axis of 150,000 single-position pieces (EC-33) under the
pan target (`itWithin`, best of three); draw and hover warm at 0 bytes on a reversed collapsed axis (the opt-in tier:
run it and report).

**Whose each choice is**: the design's table; nothing configurable beyond it; plain defaults, never inferred from what
a section draws.

**Minimal code**: extend `AxisTransform`, the path code and the hit records; reuse the codon map's runs,
`projectBoundary`, the closure, the hit's `extent` and `hull`; no private copies; sizes reported.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `pnpm test:heavy packages/seq/` (never while another heavy run is going; if one is, wait)
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`, and `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`, `pnpm format:check`
- `pnpm -C packages/seq/core size` (report the numbers)

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report);
- anything in `.seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
design, or a target cannot be met without changing a design (report the measurement and the options), or a permission
is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the sizes and the performance numbers; the
hot-file changes you need; the changelog entry's text; config items with their defaults; the matrix cells filled;
drift (what, where, options, your recommendation); followups. You cannot see a browser: name the stories to look at,
and what to look for in each (light and dark).
