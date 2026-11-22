# Lane BS: the builder showcased

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/bs`. Model:
fast (Sonnet; storybook only, well scoped). Beside VC and RG; touches no file theirs do.

**Why**: the builder (`createSeqSpec` in `@molstar/seqx-core/lang`, `lang/builder.ts`, `theme-builder.ts`,
`genome-builder.ts`; core `_spec` "Language builder"; format.md: "the builder is the main way to write a spec") is how
a host writes a view, yet no story shows it: the examples use it, the reader never sees the code. The author
(2026-11-22): "keep it in core for now, but add a storybook example to features section that showcases it".

**Task**: a story group `Features/Builder` in `packages/seq/storybook/src/_stories/features/builder/`, each story **the
builder code shown beside the view it builds** (the code a readable panel next to the canvas, copied from the same
module the story runs, so they cannot drift: import the module's source with Vite's `?raw`, or another way that keeps
one copy; say which). Short, each building one idea, a sentence of text each:

- **a first view**: a sequence from letters, one track, one layer (the smallest spec that draws);
- **handles and references**: a node's address returned and passed to a reference param (a layer over a source, a
  numbering in a sequence's body); a comment showing the line that is a type error (kept in a `// @ts-expect-error`
  so the typecheck proves it);
- **layers and themes**: `theme(...)` and `measure(...)` from `theme-builder`; a categorical and a values layer;
- **a genome view**: an assembly with the genome builder (`genome-builder.ts`), a locus;
- **check()**: a spec with a mistake the types cannot catch (an address written out that names nothing), `check()`'s
  issues listed beside the view's placeholder;
- **build() as data**: the normalized tree printed (JSON, and the S-expression printing if exported) beside the view.

Use existing fixtures (`@molstar/seqx-examples`); no network. Light and dark. Each story's headless twin (as the other
features' `_test` twins do: composed as the story composes it, on a fake clock). A row per story in
`.seq/design/features.md` and the Storybook Features page (`design/features.mdx`) are hot: list the exact rows.

**Out of scope**: no change under `packages/seq/{core,viewer,examples}` (an export a story needs: list it and stop to
ask); no builder package; no new builder API.

**Where**: a worktree on `seq/bs`. Headers on new and changed files (`@author David Sehnal <david.sehnal@gmail.com>`,
2026). Commit at each green step. Hot files never staged (`package.json`s, lockfile, `biome.json`, `.seq/`,
`.v6/changelogs/`, `storybook/packages.ts`): list exact edits.

**Checks**: `pnpm test:unit packages/seq/storybook` (under 2 s for the new folder), `pnpm -C packages/seq/storybook
typecheck`, `pnpm lint`, `pnpm format:check`, `node .seq/_checks/file-size.ts`, `node .seq/_checks/features-planned.ts
--plan .seq/execution/fifth-slice/plan.md`.

**Report**: under 250 words: commits; story ids; how the code is kept one copy; checks with numbers; hot-file edits;
the `.v6/changelogs/` text.
