# Lane C0: the second slice's carry-overs, and the bookmark and highlight commands

**Goal** (plan, stage 0; designs agreed 2026-11-18):
1. **`projectPoint`** exported from `@molstar/seqx-core/resolve` (through its existing export, no new package entry),
   and the viewer's private copy (`pathPoint` in `models/view.ts`) removed; its callers use the export.
2. **A link followed before the spec is loaded is kept waiting**, not refused: the address model treats "no spec yet"
   as waiting (as it already treats data not loaded), and a focus kept until the view is placed records its address
   when the kept move lands (today `focusAddress()` reads the shown range). The Features links story then follows its
   link through the viewer, not by holding it itself.
3. **The Features links story's single message**: one status message per followed link (not one per step of its
   resolution).
4. **`track.add`** as a command (one patch, one undo step), replacing the example's plain patch.
5. **`view.back` and `view.forward` in the default toolbar**, disabled when their history side is empty (`view.history`
   is reactive), with their keys shown in the tooltip.
6. **The selection's bookmark command** (`selection.bookmark`): adds a ref for the current selection (one patch, one
   undo step, EC-162; `addressing.md`, "Bookmarks and location sets are refs"), its key in `view.keys` (pick one with
   its own `e.key`, checked by the bindings check, e.g. Cmd/Ctrl+D), its label defaulting to the selection's printed
   address; **`refs.remove`** (a ref removed by one patch, one undo step), and a **`label` param** on a bookmark's ref,
   set by patch (the matrix's "its label and color: planned": the label here; the color is TH's).
7. **`highlight.clear`**: clears one highlight channel or all (`address.highlight`'s inverse; `themes.md`, the legend's
   highlight rows), a command with no default key.
8. **Drawing checked against the data** in the existing p53 example views: every item in view has its marks and every
   drawn mark is hit, at the fit view and across a tile edge, in a quick-lane test per view (the WRAP53 tile-edge case
   is the model).

Update `ux-foundations.md`-mirrored cells only through your report (the orchestrator owns `.scaffolds/seq/`); update the viewer
and core `_spec`s (commands, Configuration) and `packages/seq/TODO.md` (remove the items built here).

Not in scope: the axis (direction, pieces), themes, the bookmark's color, navigators, anything in `packages/core`.

**Read first**: `packages/seq/AGENTS.md`; `.scaffolds/seq/design/ux-foundations.md` (section 8, the matrix rows for the
selection and a bookmark or link), `.scaffolds/seq/design/addressing.md`, `.scaffolds/seq/design/viewer-models.md` (the history rule);
the viewer and core `_spec`s; `packages/seq/viewer/src/{commands.ts,bindings.ts,models/view.ts,models/address.ts,models/selection.ts}`,
the toolbar's chrome; `packages/seq/core/src/resolve/`; `packages/seq/TODO.md`.

**Model**: fast (sonnet).

**Base**: `packages-seq` at `c5fe373cd`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/c0-carry-overs c5fe373cd`, and check `git merge-base --is-ancestor c5fe373cd HEAD`.

**Touches**: `packages/seq/viewer/src/**`, `packages/seq/core/src/{resolve,lang}/**` (the export; the `label` param),
`packages/seq/examples/src/**` (the example's `track.add`; the views' data checks), `packages/seq/storybook/**` (the
Features links story; the toolbar), the `_spec`s, `packages/seq/TODO.md`.

**Two consumers**: a bookmark made from the selection in the viewer, read back by an SVG export of the same spec (its
ref drawn as the export draws refs) and followed as a link in a second viewer on the same spec.

**Tests**: headless twins on a fake clock: a link followed before load lands when the spec arrives and records its
focus; bookmark then undo then redo; `refs.remove` with undo; `highlight.clear` of one channel and of all; the toolbar's
disabled state across back and forward. The p53 views' drawing-against-data tests.

**Performance representative**: `selection.bookmark` under 4 ms on the p53 gene view (`itWithin`, best of three).

**Whose each choice is**: the bookmark's key (`view.keys`, `ViewerConfig`), the bookmark's default label (fixed: the
printed address, because it is what a link would say), the toolbar's buttons (`SeqViewerSpec.ui`), in the `_spec`s'
Configuration tables. Nothing configurable without a second use.

**Minimal code**: reuse what exists (the address printer, the history, the commands' patch helpers); no private copies.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `pnpm test:heavy packages/seq/`
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

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past this
brief, or a target cannot be met without changing a design (report the measurement and the options), or a permission
is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the sizes; the hot-file changes you need; the
changelog entry's text; config items with their defaults; the matrix cells filled; drift (what, where, options, your
recommendation); followups. You cannot see a browser: say what should be looked at.
