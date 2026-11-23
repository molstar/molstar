# Lane CL3: the square's way home, Hi-C in bases, the next-feature keys, the zoom hint, the square's pending state

**Model**: Sonnet (a closed item list: answers agreed gesture by gesture). **Base**: `packages-seq` at the
orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl3`. Beside you: CL10 (the novice decisions:
`viewer/src/{filter-chip,header-entries,header-parts-dom,header-icons,link-structure}.ts`,
`viewer/src/behaviors/header-extras.ts`, the structure-link stories and the fake). Never edit their files.

**Why**: the close-out plan, wave 2. The fifth slice's novice passes found nothing on the page saying the square's
window drags or resizes, a double-click moving it instead of fitting the genome, and 0 doing nothing (6); Hi-C cards in
bin indices (33); no key for the next feature (9); "zoom in to read" naming no gesture (12); the links-pending square an
empty grid (35). The answers are agreed (the author, 2026-11-23: "Agree as drafted").

**Read first**:
- `.scaffolds/seq/execution/close-out/drafts/novice-questions.md`: **C** (the gesture table is the agreed design, your
  mock), **D** item 33, **E** item 9, **G** items 12 and 35 (the square's part; the ring's hole is CL2's, merged);
- `.scaffolds/seq/execution/fifth-slice/novice-findings.md`, items 6, 9, 12, 33, 35;
- `.scaffolds/seq/design/genome-coordinates.md` (the square; A navigator over several systems),
  `.scaffolds/seq/design/ux-foundations.md` (the keys; the interaction matrix), `.scaffolds/seq/design/region-reads.md`
  (the hatch);
- `packages/seq/viewer/_spec/README.md`: Input, Bindings, the square's part (the lazy square), History; the core spec's
  issues and hatch;
- the code: `viewer/src/{square-set,square-chunk,navigator,navigators,bindings,keys,input}.ts`,
  `core/src/render/{tracks-draw,tracks-types,issues}.ts` (the hatch's words), the ring's square in `ext-rondo` where the
  square is drawn for the ring (read it; edit only what the square needs).

**Items**, a commit at each:

0. **The checks first**: `]` and `[` unbound in the binding table, every preset and the keys panel; + and = on the
   current code (the novice found them doing nothing: a real bug or the tool's keys?). A contradiction with the agreed
   answer stops the lane: report it as a drift note.
1. **The square's hint (6)**: one muted line under the square, "Drag to move · drag an edge to resize · double-click
   for the whole genome"; an option on the square's part turns it off (`hint: false`). Proof: a test; the Configuration
   table (`config-table.test.ts`).
2. **The way home (6)**: a double-click anywhere on the square gives the whole genome (replacing the move), and 0 with
   the square focused too; each one history step; Back returns the window before. Proof: a test per gesture of table C.
3. **Hi-C in bases (33)**: a Hi-C cell's card in kb or Mb ("chr1 88.5–89.0 Mb × chr1 86.5–87.0 Mb", its value), no bin
   indices; the square's rulers in kb or Mb. Proof: tests on "Features › Genome Ring › The Square Beside" (a cell's
   card) and the contact-map story.
4. **The next feature (9)**: `]` and `[` move the cursor to the next and previous feature on the focused track,
   skipping empty cells; listed in the keys panel. Proof: tests on "Features › Genome formats › GFF3 › Find By Name".
5. **The zoom hint (12)**: the hatch names the gesture ("Zoom in to read: ⌘/Ctrl+scroll, or +"); a track that draws
   nothing at this zoom gets the same hatch, never an empty band. Fix + and = if item 0 found them broken. Proof: tests on
   the Links track's case.
6. **The square pending (35)**: while links load, the square shows the hatch with "Loading…"; on a failure, the issue's
   text. Proof: tests on "Links Pending" and "Links Failed" (the square's part).
7. **The designs and specs**: the answers into `genome-coordinates.md`, `ux-foundations.md` (the keys, the matrix
   cells), `region-reads.md`, their Storybook Design pages, the viewer and core specs. Proof: `designs-in-specs` green.

**How**:
- **host neutral**: the keys go through the binding table (never a global listener); the hint is the square's part,
  removable; the history is the host's stack through the viewer's history;
- **minimal code**: reuse the navigator's jump history, the hatch, the keys panel; module size (§module-size);
- **its performance representative**: nothing per frame; the square's double-click a single fit;
- **look**: an icon in an `<svg>` gets `@molstar/ui`'s icon class or a size (jsdom cannot see a 0-by-0 icon);
- **lazy parts awaited**: the square and the ring's chunk are lazy: await them before a gesture;
- the agreed gestures (table C) are the mock; stories with headless twins on a fake clock.

**Where**: a worktree on `seq/cl3`. Touches the files above and their stories, designs, specs and Design pages. Out of
scope: CL10's files (above). Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Name them in the report (the changelog line too). The language
does not change. Write only in your worktree and your own temporary space. No network.

**Checks**, green only on exit 0 with a non-zero test count:
- `pnpm test:unit packages/seq/ cli/seqx-render/`;
- `pnpm test:medium packages/seq/`;
- `pnpm -s -C packages/seq/core size` (read the rows: **the minimal viewer has about 0.2 kB left under 1,090**; keep new
  code in the lazy square chunk or the keys panel's chunk; report every row);
- `npx tsc -b`; `pnpm lint`; `pnpm format:check`;
- the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first: `git add .scaffolds packages/seq`).

Your quick-lane time per touched folder stays under 2 s (the ring's folders are over already: add nothing slow there).
Already red at the base: name it, don't fix it.

**Stop when** the checks are green, or the same check failed twice after a fix, or a step is refused by the permission
system (report the action and why; don't work around it).

**Report** (under 400 words): commits; each item with its proof; each check with its exit code and test count; changed
expectations and why; what to look at (stories, light and dark, the gestures of table C); drift, gaps and followups; the
changelog text.
