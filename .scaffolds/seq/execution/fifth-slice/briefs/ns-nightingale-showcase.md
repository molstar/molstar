# Lane NS: a Nightingale showcase in the storybook

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/ns`. Model:
fast. Storybook only: no file under `packages/seq/{core,viewer,examples}` changes, except an export a story needs
(listed, and asked first).

**Why**: the resolved view is meant to be read by other renderers (composition.md: "other viewers (Nightingale,
JBrowse)"). NG proved it with a test-only plain-SVG consumer (EC-183). The author wants the principle shown with a real
one (2026-11-21: "yes, we can do a dev dependency version"). The fourth slice's "Nightingale dropped" stands for the
packages.

**Task**:

- add Nightingale's web components (`@nightingale-elements/*`: the sequence, track and a line or variation track, the
  smallest set that shows the three) as **dev dependencies of the Seq\* storybook package only** (hot: its
  `package.json` and the lockfile: give the exact lines and versions; the orchestrator installs). Pick versions at least
  two weeks old, check each licence (MIT expected) and say it. Biome's import rules: Nightingale imported only under the
  storybook's `_stories/`;
- a story group `Seq next / Showcases / Nightingale`, each story with a short text saying what is read from where:
  - **a sequence and its features**: p53's letters and its domains ranges read from Seq*'s resolved view and the store,
    fed to `nightingale-sequence` and `nightingale-track`, with Seq*'s own view of the same spec beside it;
  - **values**: a values layer (B-factors or a conservation track) fed to Nightingale's line track;
  - **the link**: Seq*'s hover and selection driving Nightingale's highlight (one way is enough), to show the
    resolved view's addresses carry over;
  - the NG test consumer's SVG beside them, as the reference;
- an adapter module in the story folder (`nightingale-adapter.ts`): the resolved view and store in, Nightingale's data
  shapes out, typed and small; its quick-lane test on a fixture (no DOM needed for the mapping);
- light and dark (Nightingale's colours set from Seq*'s theme tokens where it allows);
- no network: fixtures only.

**Where**: a worktree on `seq/ns`. Headers on new files (`@author David Sehnal <david.sehnal@gmail.com>`, 2026). Hot
files never staged: the storybook `package.json`, the lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`: list exact
edits. If installing in the worktree is needed to build, do it there with the exact versions you list and say so.

**Checks**: `pnpm test:unit packages/seq/storybook`, `pnpm -C packages/seq/storybook typecheck`, `pnpm lint`,
`pnpm format:check`, `pnpm -s -C packages/seq/core size` (unchanged: Nightingale in no Seq* bundle), the storybook
perf tier unchanged.

**Report**: under 250 words: commits; the dependency lines, versions and licences; story ids; checks; the
`.v6/changelogs/` text.
