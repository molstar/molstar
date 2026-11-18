# Lane N: editing in the viewer

**Goal**: build editing as `.seq/design/changes.md`, `viewer-models.md` and `patches.md`/`merge.md` (their consumers)
say, in `@molstar/seqx-viewer` (and `@molstar/seqx-core` where the model compiles increments) (plan.md, N):
- **`SpecModel.applyPatch`** with `SeqLanguage`, its validate options, and an **incremental `checkSeqLang`** as the
  patch validator (its cost in proportion to `Changes`); then resolve and the view recompiling only what `Changes`
  reaches (composition.md); Seq's own conformance cases through `checkPatches(apply, cases)`;
- **commits and recent roots**: `SpecCommit` as the port's checkpoint, `revisions.keepRecent` (default 64), the port's
  optional `root(revision)`;
- **history through Seq\*'s history port**: the port in `seqx-viewer`'s `ports`, `SpecModel` free of
  `@molstar/ui/history`, the `@molstar/ui/history` adapter in the default composition only. One step a batch, an agent's
  turn as one step, an outside change reported to the host. Tested with two participants on one host stack (a Seq\*
  spec and a stand-in structure state, undone in order), and with no history at all;
- **the rebase test**: a stale agent patch refused, then rebased by the host's opt-in through `merge3` with the base
  root from the recent roots;
- **previews as patches**, applied again when the root changes, dropped with an issue when they no longer apply (a
  preview surviving an agent's patch mid-drag), ending in one patch or nothing (Esc cancels);
- **story steps as stored patches** (`patchLanguage(SeqLanguage)`): entered and left (navigation, not history), shown
  over the base and the preview, edits while shown patching the base, masked edits reported, a step that no longer
  applies shown as such;
- **every command declares its change kind**, with the tests changes.md asks for: each command's kind matches what it
  changes, undo touches only spec patches, an arrival never changes the spec's hash.

Lane O's addressing is merged before you start: bookmarks are refs added by one patch, and its commands take
addresses. Build on them; a bookmark added or removed goes through `SpecModel.applyPatch` and history like any edit.

Not in scope: selection's behavior and drawing (lane P), the security bounds (lane Q), merging written patches (the
reserved `move`), an undo UI beyond what a story needs.

**Agreed surface**: nothing new in `@molstar/core`; `applyPatch`, `diff`, `merge3` and `patchLanguage` as they are.
Anything past them stops the lane.

**Read first**: `.seq/design/changes.md` (all of it), `viewer-models.md`, `patches.md` (Applying, the inverse and the
allowances, stored patches, previews), `merge.md` (Rebase), `composition.md`, `ports.md`, `addressing.md`; the viewer
and core `_spec`s; `packages/seq/viewer/src/{models,ports,commands.ts,viewer.ts}`, `packages/seq/core/src/{lang,view,resolve}`;
`packages/core/src/repr/_spec/README.md` (Patches, Merge); `packages/seq/AGENTS.md`.

**Model**: strong (opus).

**Base**: `packages-seq` at `61b7d3549` (lane O merged). Branch
`seq/n-editing`; commit on it only, at each small green step.

**Touches**: `packages/seq/viewer/src/**`, `packages/seq/core/src/{lang,view,resolve}/**` where increments need it,
`packages/seq/examples/src/**`, `packages/seq/storybook/**` stories under `Seq* (next)/Features/Editing` and their
headless twins; the `_spec`s. Not `packages/core`.

**What is drawn**: an editing story (hide and reorder a track, change a theme, add a bookmark; undo and redo through a
stand-in host stack; a preview during a drag that ends in one patch or Esc; the printed spec updating beside it) and a
story-steps story (enter and leave steps, edit the base while a step is shown), each with its headless twin composed
as the story composes it, on a fake clock.

**Two consumers**: an agent and the user on one history (an agent's turn as one step between the user's edits, undone
in order), and a Seq\* spec with a stand-in structure state on one host stack.

**Performance representative**: an edit (patch, resolve, redraw) on the p53 spec under 16 ms, its undo the same: a quick
case, and the heavy bound in `itWithin` (several edits after warm-ups, timed together).

**Configurable over fixed**: `revisions.keepRecent` (64), whether a stale patch is rebased (off: the host opts in),
how previews that no longer apply are reported; in the `_spec`s' Configuration tables.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `pnpm test:heavy packages/seq/`
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`, and `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`, `pnpm format:check`

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report: exports, `_test` nulls, dependencies);
- anything in `.seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
agreed surface, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the hot-file changes you need; the changelog
entry's text; config items with their defaults (and which you were unsure of); drift (what, where, options, your
recommendation); followups.
