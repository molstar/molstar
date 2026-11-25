# Lane CL10: the stage 1 novice decisions F3, F9, F11, F14

**Model**: Sonnet (a closed item list: four decisions already made, each a test). **Base**: `packages-seq` at the
orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl10`. Beside you: CL2 (the ring,
`packages/seq/ext-rondo/`, the ring's stories, some viewer header files: `header-entries.ts`, `header-extras.ts`,
`index.ts`, `viewer.ts`, `announcer.ts`) and CL4 (formats: `core/src/runtimes`, `render/row-names.ts`, the matrix,
`card-providers.ts`). Never edit their files; if a decision needs one of them, stop and report it.

**Why**: the close-out plan, wave 2. The author decided these at the fifth slice's stage 1 gate (2026-11-21: "Yes, all
five"); F25 was built then; these four were carried to this plan (the author, 2026-11-23: "Yes, one lane").

**Read first**:
- `.scaffolds/seq/execution/fifth-slice/fixes.md`: "The novice decisions, to build" (line 26, the four decisions' exact
  words) and "No track drawn as an empty well, and a sweep that holds it" (the empty-well rule F11 goes with);
- `packages/seq/viewer/_spec/README.md`: Behaviors (Headers: Entries; the "n filtered" chip), Header menu and Header
  parts, Structure link; Hover card;
- the code: `viewer/src/filter-chip.ts`, `viewer/src/header-parts-dom.ts`, `viewer/src/header-icons.ts`,
  `viewer/src/link-structure.ts`; `storybook/src/_stories/support/fake-structure/{fake-viewer,fake-view}.ts`;
  `storybook/src/_stories/features/structure-link/`.

**Items**, a commit at each:

1. **F3, "Within 5 Å" on the fake**: on a selection that is already the chain, it selects the residues of the other
   chains within 5 Å (never the selection again); its card reads "n residues within 5 Å of chain A". The fake is the
   Storybook's (no library change). Proof: a test in the structure-link stories' tests, the parent failing.
2. **F9, the colour source pressed**: the track the structure is coloured by shows its icon pressed (one at a time;
   the icon fixed in place, not moving with hover; its tooltip naming it). Proof: tests on the structure-link story
   (pressed on one track, moved when the source changes, the tooltip's text) and their twin.
3. **F11, an empty section's issue**: an empty section says its issue in its own row, with the empty-well rule (a
   one-line row naming it, never a blank well). Proof: a test on a section with no data, canvas and SVG.
4. **F14, the "n filtered" card**: it names the filter ("1 filtered: Species = Human") and says "click to clear".
   Proof: a test on the chip's card text, with two filters too.
5. **Specs**: the viewer spec's sections above updated in the commits that build them. Proof: `designs-in-specs`
   green.

**How**:
- **host neutral**: the pressed icon is the header's own part; the colour source is the structure link's state, never
  the host's;
- **minimal code**: reuse the header parts' icon slot and `@molstar/ui`'s pressed state (`aria-pressed`); module size
  (§module-size);
- **its performance representative**: nothing per frame; the header parts build when the view settles;
- **look**: an icon drawn in a `<svg>` gets `@molstar/ui`'s icon class or a size, or it lays out at 0 by 0 (CL1's
  look found exactly that; jsdom cannot);
- stories with headless twins on a fake clock.

**Where**: a worktree on `seq/cl10`. Touches `packages/seq/viewer/src/{filter-chip,header-parts-dom,header-icons,link-structure}.ts`
and the section-issue rendering in core, the structure-link stories and the fake, the viewer spec. Out of scope: CL2's
and CL4's files (above). Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s, `tsconfig*.json`,
the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Name them in the report (the changelog line too). Write only in
your worktree and your own temporary space. No network.

**Checks**, green only on exit 0 with a non-zero test count:
- `pnpm test:unit packages/seq/ cli/seqx-render/`;
- `pnpm test:medium packages/seq/`;
- `pnpm -s -C packages/seq/core size` (read the rows; the minimal viewer has under 2 kB left under 1,090);
- `npx tsc -b`; `pnpm lint`; `pnpm format:check`;
- the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first: `git add .scaffolds packages/seq`).

Your quick-lane time per touched folder stays under 2 s. Already red at the base: name it, don't fix it.

**Stop when** the checks are green, or the same check failed twice after a fix, or a step is refused by the permission
system (report the action and why; don't work around it).

**Report** (under 400 words): commits; each item with its proof; each check with its exit code and test count; changed
expectations and why; what to look at (stories, light and dark); drift, gaps and followups; the changelog text.
