# Lane S: the p53 example stories and the pages (stage 4, S1 and S2)

**Goal** (plan.md, S1 and S2): the slice's work shown as examples a reader uses, built on what is merged (the
Features stories under `Seq* (next)/Features/{Editing,Addressing,Selection}` show the pieces; these compose them on p53):
- **S1**, stories under `Seq* (next)/Examples/p53/` (`packages/seq/storybook/src/_stories/examples/p53/`), each with
  its headless twin composing it as the story does, on a fake clock, and an entry in `index.stories.ts`:
  - **Edit and undo**: on the tracks view: add a track and hide one, set a theme, undo and redo through a stand-in
    host stack, a preview by dragging a track's height (Esc cancels, release commits one step); the printed spec beside;
  - **Bookmarks and links**: on the gene's region (gene by region): a ref on a transcript that arrives later (waits,
    then resolves), a link in a short text that moves the view, "bookmark this view";
  - **Selection**: on the protein on its gene: a drag-select on the gene, projected onto the protein, in canvas and SVG.
  Use the viewer's public API and default composition as an application would; no private hooks. Each story's text
  says in one or two sentences what to try.
- **S2**: the Features page (`design/features.mdx`) and the Examples page (`design/examples.mdx`) list what is built
  now: editing (patches, undo through the host, previews, rebase), addressing (refs, links, bookmarks), selection
  (gestures, cursor, projection, column/tracks), security bounds; each row linking its story. Keep each page's style.

Not in scope: new library features; fixes beyond a story's own composition (report a bug you find instead, with a
failing test if you can).

**Read first**: `packages/seq/AGENTS.md`; the viewer `_spec`; the Features stories and their `compose.ts` under
`packages/seq/storybook/src/_stories/features/`; `examples/p53/*`; `.scaffolds/seq/design/examples.md`, `examples-and-stories.md`.

**Model**: fast (sonnet).

**Base**: `packages-seq` at the commit the orchestrator names. Your worktree may be cut from master: with a clean tree,
first `git checkout -B seq/s-stories <that commit>`.

**Touches**: `packages/seq/storybook/src/_stories/{examples,design}/**`, `packages/seq/examples/src/**` if a spec is
needed.

**Two consumers**: each story's twin composes the viewer and an SVG export (or two views) as the story does.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/storybook typecheck`
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
