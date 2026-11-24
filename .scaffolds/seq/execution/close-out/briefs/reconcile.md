# Lane RC: Finalize's reconcile (the records, not the work folder)

**Model**: Sonnet (records against the code, decided). **Base**: `packages-seq` at the orchestrator's head
(`git reset --hard packages-seq` first). Branch `seq/rc`. No lane beside you.

**Why**: `loop.md`, Finalize, step 1: the designs, the specs, `status.md`, the features rows, the edge cases' outcomes
and version 1's rows agree with what the close-out built. The work folder is the orchestrator's (never edited here).

**Read first**: `.scaffolds/seq/execution/close-out/plan.md`, `checklist.md`, `findings.md` (every lane's row: what
it built, its commit), `drafts/version-1-rows.md` (each row's decision: build, defer with its trigger, cut, and the
author's "Agree all" in `gaps.md`); `.scaffolds/seq/status.md`; `.scaffolds/seq/design/features.md`;
`.scaffolds/seq/design/edge-cases.md` ("Outcomes of …" sections, their shape) and
`.scaffolds/seq/_checks/ec-ids.ts` (what `--worklist` reads).

**Items**, a commit at each:

1. **The features rows**: every row version 1's decisions named moves to what the close-out did, with the evidence (a
   story, a test file, a design): built rows `done` (or `building` with the part not built named); deferred rows
   `deferred` with their trigger as decided; cut rows say so. Rows the close-out's fix lanes and the novice pass's
   gaps added to (search, track groups, values marks, the legend's layers menu, `height: 'content'`, the ring's
   gestures) say it. Proof: the features check in `features-planned` and `designs-in-specs` green.
2. **The edge cases' outcomes**: `node --experimental-strip-types .scaffolds/seq/_checks/ec-ids.ts --worklist` lists
   the built rows' items no test names (56 lines at the orchestrator's head). For each: a test that names it already
   exists under another wording (cite it, `EC-<n>` in the test's name or comment), or an outcome row in
   `edge-cases.md`'s close-out section ("Outcomes of the close-out": kept as designed, changed and how, or tested by
   which file), never a test written only to cite a number. Proof: `--worklist` empty, or every item listed with an
   outcome; `ec-ids` green.
3. **`status.md`**: the current state after the close-out: the language at `0.11.0` (groups, reference lines,
   `height: 'content'`), what was built (search, track groups, values marks, the legend's layers menu with Show and
   Hide all, the ring's double-click and gestures, `track.height`), what is open (the TODO's new lines, the carried
   split and merge bound, the old Seq\* viewer's label churn), what comes next (the next plan, unwritten: say so).
   Current state only. Proof: `links` green.
4. **The designs and specs**: read each design the close-out touched (`search.md`, `track-groups.md`,
   `values-marks.md`, `themes.md`, `track-header-parts.md`, `ring-extension.md`, `layouts.md`, `across-sections.md`,
   `tooltips.md`, `performance.md`) against its Storybook page and the packages' `_spec`s; fix a disagreement in the
   design or page (never the code), and list each. Proof: `designs-in-specs` green; the list.

**Where**: a worktree on `seq/rc`. Headers per the root `AGENTS.md` (only for source files you touch; Markdown has
none). Never stage the work folder, `packages/seq/TODO.md`, the hot files (`package.json`s, `tsconfig*.json`, the
lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/_measure/sizes.mjs`). No code
changes except a test's name or comment citing an edge case. Never a bare `git stash`. No network.

**Checks**, green only on exit 0: `pnpm format:check`; `pnpm test:unit packages/seq/` (if a test's name changed); the
scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first), `ec-ids` and `features-planned --plan
.scaffolds/seq/execution/close-out/plan.md` among them.

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 350 words): each item's changes and its proof; the rows moved (counts: done, deferred, cut); the
edge cases (named by a test, given an outcome); the design and spec disagreements found and fixed; the checks.
