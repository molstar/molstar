# Lane CL2: the ring's stem menu (route, partners), its cards, links pending

**Model**: strong (a design's core: the stem menu's structure and the ring's input). **Base**: `packages-seq` at the
orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl2`. Beside you: CL14 (scaffold checks) and
CL1 (the legend, `packages/seq/{core,viewer}`'s legend code). Never edit `viewer/src/legends.ts`, `mount-legend.ts` or
the themes design; CL1 owns them.

**Why**: the close-out plan, wave 1. The fifth slice's novice passes found the stem menu unreachable (28), its partners
repeated with no word on what separates them (8, 26), the stem card's first row misread (21), a pair's card naming one
end (25), and the links-pending ring saying nothing (35). The answers are agreed (the author, 2026-11-23: "Agree as
drafted").

**Read first**:
- `.scaffolds/seq/execution/close-out/drafts/novice-questions.md`: **B** (the gesture table is the agreed design, your
  mock), **D** (21, 25) and **G** (35, the ring's part);
- `.scaffolds/seq/execution/fifth-slice/novice-findings.md`, items 8, 21, 25, 26, 28, 29, 35 (what the novice saw);
- `.scaffolds/seq/design/ring-extension.md`: Stems; Ribbons; Hover, selection, emphasis, the keyboard; The hole; The
  interaction matrix's cells; Configuration;
- `.scaffolds/seq/design/tooltips.md` (the card composed by providers);
- `packages/seq/ext-rondo/_spec/README.md` (Stems, Ribbons, Hits, emphasis, words; Configuration) and the viewer spec's
  Hover card and Header menu sections;
- the code: `packages/seq/ext-rondo/src/{menu,partners,card,ring,viewer,words}.ts`.

**Items**, a commit at each:

0. **The checks first**: find what a "stem piece" is in `RondoLay` (the stems a chromosome splits into) and whether
   each piece is one contiguous run, so a band range can name it; confirm item 29 (Shift+F10 after a click opens
   nothing) on the current code. Anything that contradicts the agreed answer stops the lane: report it as a drift
   note, don't redesign.
1. **The route (28, 29)**: a right-click on a stem opens the stem menu at the pointer; over the ring's empty space the
   browser's own menu still opens. Shift+F10 and the Menu key on the focused stem open it at the stem. Proof: a test per
   gesture of table B, by name.
2. **Partners by chromosome (8, 26)**: one row per partner chromosome with its total links, by count; a chromosome
   reached by more than one ribbon gets a submenu of its pieces, each named by band range and strand ("q34.1 · 375
   links", "q21–q22 · 24 links, reversed"). Proof: tests on the menu's rows for the "Kinds Of Links" fixture, chr9's
   three rows and chr22's two.
3. **Links within**: the stem's own chromosome reads "Within chr22 · n links", listed last, never as a partner. Proof:
   a test.
4. **Only the stem's actions**: Move up, Reverse and Hide track leave the stem menu (they stay in the track's header
   menu). Proof: a test on the menu's items.
5. **The stem card (21)**: the first row names what it counts ("This stem: 24 links to chr9"); the per-layer totals
   under a heading "chr22, all links". Proof: a test on the card's rows.
6. **A pair's card (25)**: both ends, each with its position and band; the bins row "n bins: zoom in to see them".
   Proof: a test on "A Pair"'s ribbon card.
7. **Links pending on the ring (35)**: while the links load, the hole says "Loading links…" (muted, centred), announced
   once, politely; on a failure the issue's text replaces it. The square's pending state is CL3's: leave it. Proof:
   tests on "Links Pending" and "Links Failed" and their twins.
8. **The designs and specs**: the answers into `ring-extension.md` and `tooltips.md`, their Storybook Design pages,
   and the ext-rondo spec (and the viewer spec where the card's providers say so), in the commits that build them; the
   interaction matrix's cells for the stem menu filled. Proof: `designs-in-specs` green; the pages updated.

**How**:
- **host neutral**: the right-click is the ring's own over a stem only; the browser's menu elsewhere. The keymap stays
  the viewer's binding table (no global listener); the host's chrome is untouched;
- **two consumers**: the partner list is read by the menu and the partner keys (`ThePartnerKeys`): test both live in
  one test;
- **minimal code**: reuse `partnersOf` and the header menu's submenu; no second menu implementation; module size
  (§module-size): a file past 800 lines is a question in the report;
- **its performance representative**: building the menu for the 3,000-contig ring ("Crowded") stays out of the
  frame; a quick test that its partners list for a stem is bounded;
- **lazy parts awaited**: the ring's chunk and the square are lazy; await them before a gesture;
- **the agreed gestures** (table B) are the mock; stories with headless twins on a fake clock (the loop, Each step 8).

**Where**: a worktree on `seq/cl2`. Touches `packages/seq/ext-rondo/`, the viewer's hover-card providers only where the
card needs it, `packages/seq/storybook/src/_stories/` (the ring's stories, twins, Design pages),
`.scaffolds/seq/design/{ring-extension,tooltips}.md`. Out of scope: the legend (CL1), the square's code (CL3). Headers
per the root `AGENTS.md`. Hot files never staged: `package.json`s, `tsconfig*.json`, the lockfile, `biome.json`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`, the work folder.
Name them in the report (the changelog line too). The language does not change. Write only in your worktree and your own
temporary space. No network.

**Checks**, green only on exit 0 with a non-zero test count:
- `pnpm test:unit packages/seq/ cli/seqx-render/`;
- `pnpm test:medium packages/seq/`;
- `pnpm -s -C packages/seq/core size` (read the rows; the ring's row is under 45 kB);
- `npx tsc -b`; `pnpm lint`; `pnpm format:check`;
- the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first: `git add .scaffolds packages/seq`).

Your quick-lane time per touched folder stays under 2 s. Already red at the base: name it, don't fix it.

**Stop when** the checks are green, or the same check failed twice after a fix, or a step is refused by the permission
system (report the action and why; don't work around it).

**Report** (under 500 words): commits; each item with its proof; each check with its exit code and test count; changed
expectations and why; what to look at (stories, light and dark, the gestures of table B); drift, gaps and followups; the
changelog text.
