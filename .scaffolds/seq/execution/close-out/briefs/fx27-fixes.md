# Lane FX27: the ring's Cmd+wheel and the 1M lookup's bound

**Model**: Sonnet (a closed item list, each decided). **Base**: `packages-seq` at the orchestrator's head
(`git reset --hard packages-seq` first). Branch `seq/fx27`. Beside you: CL5 (search: a new `core/src/search/`,
`viewer/src/models/` search and highlight code, `viewer/src/config.ts`'s keys) and CL6 (track groups: `core/src/lang/`,
`core/src/resolve/`, the tracks layout and headers, `viewer/src/models/view.ts`, the headers' behaviors). Never edit
their files; the specs (`_spec/*.md`) and `viewer/src/config.ts` are shared: add to your own sections and entries only.
**`viewer/src/models/view.ts` is CL6's**: if item 1 needs it, stop and report. FX25's parked worktree holds
`core/src/render/` ranges and lanes and the variant card: avoid those.

**Why**: the author's decisions, 2026-11-24 ("Agree all"): `drift/ring-cmd-wheel-at-fit.md`, option (c); and
`fixes.md`'s "The 1M-range lookup's ratio fails on a quiet machine too" (the TODO's rule: measure at a larger scale,
never raise to pass).

**Read first**: `.scaffolds/seq/execution/close-out/drift/ring-cmd-wheel-at-fit.md` (the evidence, option (c));
`.scaffolds/seq/design/ring-extension.md` (The ring zoom: a cut; its matrix), `layouts.md` (Rings: a ring nobody follows
is idle), `ux-foundations.md` (the wheel, the matrix row "a ring section"); the viewer spec's Input, Behaviors
(`navigate`) and History; `viewer/src/behaviors/navigate.ts` (`wheel`), `viewer/src/input.ts` (`wheel`),
`ViewModel._idleRing` (FX26: read only), `ext-rondo/src/menu.ts` (the stems' `section.include` args, `onlyStem`);
`core/src/data/_test/intervals.test.ts` and the interval lookup it times; `packages/seq/TODO.md`'s "M5's leftovers"; the
root `AGENTS.md`'s Tests (heavy lane, `itWithin`, `test-harness/src/timing.ts`).

**Items**, a commit at each:

1. **Cmd/Ctrl+wheel over an idle ring opens a chromosome** (option (c)): in over a chromosome's block, it is opened as
   a cut (`section.include` of that block, the args the stem menu's "Show only" builds), one patch, one undo step; out
   over an opened chromosome shows all again; one cut a wheel gesture (deltas past a notch, then nothing until the
   wheel rests, `coalesceWheel`'s rest); over the hole (no block under the pointer), nothing changes and the first such
   wheel of a gesture says once, in the live region, "The ring shows the whole genome: open a chromosome from its
   label's menu or the square". A ring someone follows keeps today's zoom. The ring's knowledge stays where the ring's
   code is (core's rings layout or `ext-rondo`), the viewer's input generic. Proof: tests per direction (in over a
   block, out over an opened one, over the hole, a followed ring unchanged), one undo step each, the event taken in
   each; the Genome Ring stories' headless twin.
2. **The 1M-range lookup measured by its growth**: replace "within 1.5× of 10,000" by the lookup's growth across
   10⁴, 10⁵, 10⁶ and 10⁷ ranges (log-like: each tenfold step under a ratio the design's interval index promises, the
   ratio from the measured best of three, with its margin said), under `itWithin` in the heavy lane; the 10⁷ fixture
   built once a file. If the measured growth is not log-like, stop and report the numbers: that is a finding, not a
   bound to tune. Proof: the new heavy case green alone three times in a row (say the load each time; under load 3
   the numbers are indicative).
3. **The records**: the viewer spec (Input, Behaviors) and `ring-extension.md`'s matrix cell, `ux-foundations.md`'s
   row "a ring section" (zoom: the cut) with their Storybook Design pages; `testing.md` or `performance.md` where the
   lookup's bound is stated. Proof: `designs-in-specs` green.

**How**: host neutral (the wheel taken only where something happens or is said; Ctrl+wheel over the hole on the page
still taken, said once); nothing per frame; module size (§module-size: `input.ts` is 1,421 lines, add nothing there
that can live in the behavior).

**Where**: a worktree on `seq/fx27`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm test:heavy packages/seq/core/src/data/_test/intervals.test.ts` alone;
`pnpm -s -C packages/seq/core size` (every row: headless 978.7 of 980, default 1,557.3 of 1,560, minimal 1,096.6 of
1,100, ring 45.0 of 48); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the scaffold checks in
`.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report the
action and why; don't work around it, not even by another tool or the Edit tool).

**Report** (under 350 words): commits; each item with its proof; item 2's numbers with the load; the checks with exit
codes and counts; the size rows; what to look at (the Genome At The Fit, Cmd+wheel in, out and over the hole, light and
dark); drift, gaps and followups; the changelog text.
