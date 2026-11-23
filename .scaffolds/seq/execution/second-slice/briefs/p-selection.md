# Lane P: selection, highlights and the cursor

**Goal** (plan.md, P), on lane O's minimal `SelectionModel` and `HighlightModel` (`packages/seq/viewer/src/models/{selection,address}.ts`):
- `SelectionModel`: element keys and a Location, the anchor, the cursor; the selection holding its nodes (dynamic nodes
  stay while selected); selection as addresses (lane O's `selectionAddress`, `followMoves` kept working);
- `HighlightModel`: named channels (lane O's `link` and `link-hover` among them) drawn in the `highlight` pass;
- the default preset's click, Shift-click and drag-select bindings and the keyboard cursor; the grammar's toggle
  (ux-foundations.md section 8);
- `Emphasis.projected`: the hover and the selection projected onto linked sections and drawn (the first slice's open
  item), the same look lighter;
- the selection tint drawn over the data (the first slice's decision), its place a config item;
- the Links story's highlights drawn (it says now they are listed, not drawn: update its text).

Not in scope: the security bounds (lane Q), story steps, anything in `packages/core`. Lane Q edits `render/svg.ts` for
escaping at the same time: leave that file's text output alone.

**Read first**: `.scaffolds/seq/design/ux-foundations.md` (sections 6 to 8), `viewer-models.md`, `addressing.md`,
`locations-and-queries.md`, `layout-and-draw-lists.md` (passes); the viewer and core `_spec`s;
`packages/seq/viewer/src/{models,behaviors,bindings.ts,commands.ts}`, `packages/seq/core/src/render`; `packages/seq/AGENTS.md`.

**Model**: strong (opus).

**Base**: `packages-seq` at `3e6f9267d`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/p-selection 3e6f9267d`, and check `git merge-base --is-ancestor 3e6f9267d HEAD`.

**Touches**: `packages/seq/viewer/src/**`, `packages/seq/core/src/{render,query,view}/**` where emphasis needs it,
`packages/seq/storybook/**` (a selection story under `Seq* (next)/Features/Selection` on both p53 views, the tracks
and the protein on its gene, with its headless twin on a fake clock), the `_spec`s.

**What is drawn**: selection, hover projected onto the linked section, highlights on their channels, the cursor; in
canvas and SVG alike (the backends' tests), light and dark.

**Two consumers**: the selection made in the viewer drawn in an SVG export of the same spec, and projected onto the
gene section in both.

**Performance representative**: a drag-select step and a cursor step under 4 ms (`itWithin`, several after warm-ups);
hover and selection redraws allocation-free warm (the opt-in tier, `pnpm -C packages/seq/core test:perf`).

**Configurable over fixed**: the tint's place and opacity, the projected look's strength, the cursor's step keys, in
the `_spec`s' Configuration tables.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `pnpm test:heavy packages/seq/`
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`, and `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`, `pnpm format:check`

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.scaffolds/seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.scaffolds/seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults (and which you were unsure of); drift (what, where, options, your
recommendation); followups.
