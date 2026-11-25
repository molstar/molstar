# Lane SO: Features apart from Examples

**Goal**: the rule now in [examples-and-stories](../../../design/examples-and-stories.md) ("Features and examples are
apart", the author 2026-11-19: "we should distinguish functionality/feature explanation from examples"): a story that
explains how something works is `Seq* (next)/Features/<Feature>/<Story>`, on whichever dataset shows it best; a
dataset's figure is `Seq* (next)/Examples/<Dataset>/<View>`. A dataset's index story links the features shown on its
data.

Today `packages/seq/storybook/src/_stories/examples/p53/` holds stories of both kinds. Move to Features, at least:
Themes (with Editing Themes and the Schemes Gallery), Edit and undo, Bookmarks and links, Selection, Story; and any
other story whose point is a function rather than p53's figure (judge each: Text, Genomic follower, Reversed gene,
One axis, the strips; say why for each kept or moved). `features/` already has addressing, editing, layers,
navigators, selection: merge with them where they overlap (one story per thing, no duplicates), keeping the better of
two. Keep each story's headless twin with it (the twins in `examples/_test/` follow their stories; their tests keep
passing, with the same coverage). Update every link: the Design pages (`_stories/design/*.mdx` `?path=/story/...`
links), the p53 index story, `example-root.ts`, any `describe` names that cite a story path. The examples Design page
lists the Features beside the Examples.

**Not in scope**: anything else; no behaviour change. Lane FS's design is under review, no code lane runs beside you.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/so-stories-split
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/` (the same test count as the base, or more; say if fewer and why);
`pnpm test:medium packages/seq/`; `pnpm -C packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; and a
link check: every `?path=/story/seq-next-...` in the repository's `packages/seq` resolves to a story id that exists
(write it as a test in the storybook package's `_test` if none exists, so a moved story can't leave a dead link).

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a moved file keeps its header.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, every `package.json` and `tsconfig*.json`,
anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix.

**Report**: commits; the new tree (Features and Examples, story by story, moved or kept and why); checks with exit codes
and test counts; links updated.
