# Lane CL1: a layer's legend (hiding it, seeing it hidden), synteny's legend, legends' domains

**Model**: strong (a design's core: it changes the legend's rule in the viewer spec, and every view with a legend uses
it). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl1`.
Beside you: CL14 (scaffold checks) and CL2 (the ring, `packages/seq/ext-rondo`). Never edit `ext-rondo` files; CL2 owns
them.

**Why**: the close-out plan, wave 1. The fifth slice's novice passes found no way to hide a kind of link from its
legend (5), no way to see which layer is hidden (13), and no legend on synteny (7). The answers are agreed (the author,
2026-11-23: "Agree as drafted"). Feature row 113 "legends state their domain" is planned for version 1 (the author,
2026-11-23: "take the rest as recommended").

**Read first**:
- `.scaffolds/seq/execution/close-out/drafts/novice-questions.md`: **A** (the gesture table is the agreed design, your
  mock) and **E**, item 7;
- `.scaffolds/seq/design/themes.md`, "Legends: data in the core, a part a host places";
  `.scaffolds/seq/design/across-sections.md` (A bridge; Defaults; Configuration);
- `packages/seq/viewer/_spec/README.md`: Legends (the click rule, "Hiding has two owners", Hidden tracks), Commands,
  Bridges, Configuration;
- the code: `viewer/src/{legends,mount-legend,header-entries,commands}.ts`, `viewer/src/behaviors/theming.ts`,
  `core/src/render/bridges.ts`, the core theme's legend derivation (`core/src/theme/`).

**Items**, a commit at each:

0. **The checks first** (the answer's **Check** items): confirm `layer.hide` and `layer.show` are patches, one undo step;
   find whether the bridges' colour is a channel (if fixed in the draw, it becomes one, so its legend derives as every
   other one); find what "legends state their domain" lacks today. Anything that contradicts the agreed answer stops
   the lane: report it as a drift note, don't redesign.
1. **A layer's eye**: a legend row for one layer gets an `@molstar/ui` icon button after its title ("Hide <layer>");
   a merged legend over several layers gets one button opening a menu of its layers. It runs `layer.hide`. Proof: a
   test per gesture of table A, by name.
2. **A hidden layer's row stays**: dimmed title and swatches, the button "Show <layer>" (eye-off, `layer.show`), its
   entries with no action and a not-allowed pointer. Undo restores in one step. Proof: tests by name; the story
   "Features › Genome Ring › Kinds Of Links Contacts Hidden" (and its twin) shows the dimmed row.
3. **The host's option**: `mountLegend(…, { layers: false })` drops the eye (on by default). Proof: a test; the
   Configuration table updated (`config-table.test.ts` passes).
4. **Look 24**: the Contacts entry gets its swatch; its label stops doubling ("Contacts: by chromosome"). Proof: a test
   on the entry's swatch and label.
5. **Synteny's legend (7)**: the bridges' colour by orientation (same, inverted) as entries, and a note that a bridge's
   width is its block's length; an inversion's colour differs from a same-strand bridge's. Proof: a test on
   `viewLegends` for the synteny view; "Examples › Genomes › Views › Synteny" shows it.
6. **Legends state their domain (113)**: a continuous legend over an open domain says its domain (the settled min and
   max, as the row's "open domains at settle" says). Proof: a test; the feature row moves to built (item 8).
7. **The designs and specs**: the answers into `themes.md` (Legends) and `across-sections.md`, their Storybook Design
   pages, and the viewer spec's Legends and Bridges sections, in the commits that build them. The spec gains "a layer's
   row keeps its place when hidden". Proof: `designs-in-specs` green; the pages updated.
8. **The features rows**: 113 moved to built in `features.md` and the Storybook's `design/features.mdx`. Proof: the
   rows.

**How**:
- **host neutral**: the eye is in the legend part a host places; it never adds chrome to the view, and editing stays
  the application's (`layer.hide` is an existing session-reachable patch, not a new owner of the legend);
- **two consumers**: the eye's state is tested with the legend part and the header menu's "Hidden layers" live in one
  test (they must agree);
- **minimal code**: reuse `@molstar/ui`'s icon button and menu; no new legend model, the derivation stays in core;
  module size (§module-size): `mount-legend.ts` past 800 lines is a question in the report;
- **its performance representative**: legends derive when the view settles, never per frame; the existing "legends
  for 100 tracks under 1 ms" bound stays green;
- **lazy parts awaited**: the ring's stories have lazy parts; await them before a gesture;
- **the agreed gestures** (table A) are the mock: every row a test by name; stories with headless twins on a fake
  clock (the loop, Each step 8).

**Where**: a worktree on `seq/cl1`. Touches `packages/seq/{core,viewer}/src`, their `_spec/README.md`,
`packages/seq/storybook/src/_stories/` (stories, twins, Design pages), `.scaffolds/seq/design/{themes,across-sections,features}.md`.
Out of scope: `packages/seq/ext-rondo/` (CL2), the track header's other parts (CL10). Headers per the root `AGENTS.md`.
Hot files never staged: `package.json`s, `tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`, the work folder. Name them in the
report (the changelog line too). The language does not change: a change to it stops the lane. Write only in your
worktree and your own temporary space. No network.

**Checks**, green only on exit 0 with a non-zero test count:
- `pnpm test:unit packages/seq/ cli/seqx-render/`;
- `pnpm test:medium packages/seq/`;
- `pnpm -s -C packages/seq/core size` (read the rows);
- `npx tsc -b`; `pnpm lint`; `pnpm format:check`;
- the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first: `git add .scaffolds packages/seq`).

Your quick-lane time per touched folder stays under 2 s. Already red at the base: name it, don't fix it.

**Stop when** the checks are green, or the same check failed twice after a fix, or a step is refused by the permission
system (report the action and why; don't work around it).

**Report** (under 500 words): commits; each item with its proof; each check with its exit code and test count; changed
expectations and why; what to look at (stories, light and dark, the gestures of table A); drift, gaps and followups; the
changelog text.
