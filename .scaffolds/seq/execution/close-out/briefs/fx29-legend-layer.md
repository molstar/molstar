# Lane FX29: a legend row that reads as one category and hides the whole track

**Model**: Sonnet. **Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch
`seq/fx29`. Beside you: FX28 (a bridge's ribbons on the canvas: `core/src/view/project.ts`, the bridge's lay and draw
in `core/src/render/bridge*.ts`) and CL8 waiting to merge (the values layer, `values-*.ts`). Never edit their files; one
you need: stop and report.

**Why**: the author, 2026-11-24, two screenshots of Features › Structure link › Focus And Back: "there is a bug on the
legend here, showing just the ATP, but also hiding everything". The legend's row reads "Sequence ⊙ ATP" (the track's
name and one red swatch) while the Sequence track draws its letters, ATP ×2 in red and six other ligands (PYZ, 2AN,
B49, STU, LZ1, 4SP) in grey; clicking the row's eye empties the whole track (letters and every ligand gone) and leaves
the header truncated to "S…" with a stray square marker.

**Read first**: `.scaffolds/seq/execution/close-out/drafts/novice-questions.md` answer A (a layer's legend: the eye,
"Hidden, the row stays", the gesture table) as built by CL1; `.scaffolds/seq/design/themes.md` (Legends: what a
layer's row lists, rules and categories), `track-header-parts.md` (the header's label and its parts); the viewer spec's
Legends; `viewer/src/legends.ts`, `legend-eye-dom.ts`, `mount-legend.ts`, `core/src/render/legends.ts`; the story's
spec in `storybook/src/_stories/features/structure-link/`.

**Items**, a commit at each:

0. **Find what the row is**: which layer(s) the "Sequence" row stands for (the letters layer with a theme rule for
   ATP, a separate ligand layer, or a merged row), why it lists ATP alone (the other ligands and the letters' own
   colouring left out), and what its eye hides. If the designs say what such a row lists or hides and the code
   disagrees, that is the bug; if the designs do not say (a letters layer whose theme has one highlight rule over a
   base colouring), stop that part with a drift note in your report (options and a recommendation), and go on with
   the rest. Proof: the finding in the report, with the spec path.
1. **The row names what it hides**: a layer's row is titled by the layer (its label, else its kind and the track's
   name: "Sequence letters"), and its eye's label says the same ("Hide Sequence letters"), so the row never reads as
   one category's switch; a merged row's eye opens the menu of its layers (answer A). Proof: a legend test on this
   story's composition.
2. **The row lists what the layer draws**: every category its theme draws in view (the ligands in grey as their
   category, the letters' scheme as its own entry or a "+ residue colours" line, whatever the design or the drift's
   recommendation says), never only the rule's. Proof: the same test.
3. **The header after a hide**: the track's header keeps its label (no "S…" truncation from a part that appears only
   when hidden) and the stray square is explained or removed (a hidden-layer marker drawn where the design puts it).
   Proof: a header layout test on the hidden state.
4. **The class**: a sweep over every story with a legend (medium lane) asserting each row's title names its layer and
   its entries cover the layer's drawn categories; a quick representative. Proof: the sweep.

**How**: host neutral; nothing per frame; module size; reuse the legend's existing model.

**Where**: a worktree on `seq/fx29`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network. Do not run the heavy lane.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` and `pnpm -s -C packages/seq/storybook
test:perf`; `pnpm -s -C packages/seq/core size` (every row); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the
scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 350 words): item 0's finding; commits; each item with its proof (failing on its parent); the checks
with exit codes and counts; the size rows; what to look at; drift; the changelog text.
