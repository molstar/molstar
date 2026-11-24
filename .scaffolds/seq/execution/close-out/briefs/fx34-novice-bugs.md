# Lane FX34: the close-out's novice pass, its bugs

**Model**: Sonnet (bugs, each measured first). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard
packages-seq` first). Branch `seq/fx34`. No lane beside you.

**Why**: [novice-findings.md](../novice-findings.md), the rows placed to FX34 (reproduced by the orchestrator where
marked).

**Items** (measure first, then fix, each with a test failing on its parent; a commit each):

1. **C5, a partner from the stem menu is no history step** (reproduced: Features › Genome Ring › The Stem Menu, "Show
   with a partner" › chr12 › its first piece: the pair opens, Back stays `aria-disabled`). `ring-extension.md`, The stem
   menu: "Click a partner (or a piece): shows the pair, one history step". Make it one (Back returns to the ring as it
   was; one step, not two). Test headless through the menu's model.
2. **B4, A8, the stories' SVG twin misses state** (reproduced: Examples › p53 › Search, Clear: the canvas's tints gone,
   the SVG twin's kept; Features › Track groups › Opening, a fold then an undo: the twin open, the canvas folded). Find
   what the story's twin re-exports on (the composition in `storybook/src/_stories/`, the export's inputs: channels, the
   fold state) and make it follow every state the canvas draws; fix the class (every story using the twin), not two
   stories. Test: the twin's export after a channel clear and after an undo of a fold equals a fresh export.
3. **A6, a folded group's summary layer named by its raw id** (`packages/seq/TODO.md`, FX29's followup: "ranges#summary",
   "secondary#summary" in its legend row and eye; right after a fold its own name shows, then the id). Name it as the
   group's track it summarizes (its label, else its kind and the group's label), the same in the row, its eye and the
   hidden row; and find why the name changes after the fold settles. Test on Features › Track groups' composition.
4. **A4, "letters letters"** (Features › Legends › Long Legend): a track labelled "letters" titles its letters layer's
   row "letters letters". The layer's name rule (`layerName`, FX29) must not double a word the label already ends with.
   Test.
5. **B7, B10, Find By Name's lanes and its line**: Features › Genome formats › GFF3 and GTF › Find By Name: their
   Transcripts track at `height: 'content'` (FX30's param), so canvas and SVG draw the same lanes; the result block's
   first line wraps instead of being cut (GTF, dark). Proof: the stories' twins.
6. **B12, Manhattan's card**: "chr 11,984, index 11983, -log10 p 5.78": a values card shows the position (as every
   other card does), never an index. Test on the story's composition.

**Records**: the specs a rule moves (the ring's stem menu if it changes, the legend's names), `examples.md` if a story
moves; `packages/seq/TODO.md` is the orchestrator's (A6's line is removed at the merge).

**How**: nothing per frame; host neutral; module size.

**Where**: a worktree on `seq/fx34`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network. Do not run the heavy lane.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` (load under 3: `uptime` first);
`pnpm -s -C packages/seq/core size` (every row: spec 265.0 of 268, headless 995.0 of 1,000, default 1,588.1 of 1,592,
minimal 1,119.7 of 1,124); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the scaffold checks in
`.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 350 words): each item's cause; commits; each item with its proof (failing on its parent); the checks
with exit codes and counts; the size rows; the stories to look at; the changelog text.
