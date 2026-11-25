# Lane CL14: the scaffold checks (a brief's items; the heavy bounds isolated)

**Model**: Sonnet (a closed item list and named outputs: one check mode and one config change). **Base**:
`packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl14`. Beside you: CL1
(the legend) and CL2 (the ring), in other files.

**Why**: the close-out plan, wave 1 (`plan.md`, Lanes). These two are the fifth slice's fold-back B3 and B7, decided by
the author on 2026-11-23 ("Apply all ten"), and left for this plan:
- **B3**: a brief's items are checked as landed, not as reported. GC2 merged with `section.center` in its brief and
  never built it, and nothing compared the report with the brief.
- **B7**: the heavy lane's tight bounds fail under the lane's own worker contention and pass alone. Isolate them; never
  raise them.

**Read first**:
- `.scaffolds/_method/_checks/features-planned.ts` and its `--self-check`; `.scaffolds/_method/_checks/lib.ts`;
- `.scaffolds/_method/brief.md` (the Items section) and `.scaffolds/_method/findings.md` (the row);
- `.scaffolds/_method/loop.md`, Each step 5, "items landed, not reported";
- `vitest.config.ts` (the lanes, `Isolated`, the heavy mode's workers); `test-harness/src/timing.ts` (`itWithin`);
- `.scaffolds/seq/principles/README.md`, Adaptations, {#heavy-lane-alone};
- `packages/seq/TODO.md`, the line on the three heavy bounds under contention.

**Items**, a commit at each:

1. **A brief's items against its findings row**: `features-planned.ts` gains `--brief <brief.md> --findings
   <findings.md>`. It reads the brief's numbered items (`1. **<name>**`) and the lane's row in `findings.md` (by the
   brief's lane id). It fails when an item is neither ticked with a proof nor marked "not landed" in the row, and names
   each one. Define the Items cell's form (for example `1 ✓ <proof>; 2 ✗ not landed: <where it went>`) and write it into
   `findings.md`'s template and `loop.md`'s record step. Proof: the check on a planted brief and row, both ways.
2. **Its self-check**: planted cases in `--self-check`: a row ticking every item passes; an item missing, an item
   ticked without a proof, and a row missing each fail. Proof: `features-planned.ts --self-check` exit 0, the planted
   failures named in its output.
3. **The three bounds isolated**: the 1,000,000-range lookup (`core/src/data/_test/intervals.test.ts`), the pairs
   pyramid's slices (find it: a pairs or region test under `itWithin`), and the field theme's 1.1× ratio
   (`core/src/render/_test/themes.test.ts`). Each runs alone in its worker in the heavy lane: a project or a file group
   in `vitest.config.ts` that runs after the rest, one file at a time. Say why in a comment there. No bound is raised and
   no `tries` added to hide it. Proof: `pnpm test:heavy packages/seq/` green, with the three tests named in the output
   of the isolated group.
4. **Every bound under 50 ms isolated the same way** (the manifest's rule): list the heavy `itWithin` bounds under
   50 ms in `packages/seq` and put each in the isolated group, or say in the report why one needs it not. Proof: the
   list in the report and the config's group.

**How**:
- **minimal code**: reuse `lib.ts`'s argument and finish helpers; the new mode is a function beside `judge`.
- **the check stays the method's**: it reads only the method's brief and findings shapes, never a Seq\* path.
- **host neutral, lazy parts, the mock**: not applicable (no drawing).

**Where**: a worktree on `seq/cl14`. Touches:
- `.scaffolds/_method/_checks/features-planned.ts`, `.scaffolds/_method/{findings,loop,brief}.md`;
- `vitest.config.ts` (allowed for this lane: it is not on the hot list);
- the three tests' files, only to move or tag them, never their bounds.

Out of scope: any other check, any package's code. Headers: the root `AGENTS.md` rule on every source file you change
(the `.ts` files; Markdown has none). Hot files never staged: `package.json`s, `tsconfig*.json`, the lockfile,
`biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`,
the work folder (`.scaffolds/seq/execution/close-out/`). Name what you need in the report. Write only in your worktree
and your own temporary space. No network (`pnpm install --offline` if needed).

**Checks**, each its own command from the worktree's root, green only on exit 0 with a non-zero test count:
- `node .scaffolds/_method/_checks/features-planned.ts --self-check`;
- `node .scaffolds/_method/_checks/features-planned.ts --instance .scaffolds/seq --plan .scaffolds/seq/execution/close-out/plan.md`;
- `pnpm test:unit packages/seq/core`;
- `pnpm test:heavy packages/seq/`, alone (no other test run beside it);
- `pnpm lint`; `pnpm format:check`.

Already red at the base: name it, don't fix it.

**Stop when** the checks are green, or the same check failed twice after a fix, or a step is refused by the permission
system (report the action and why; don't work around it).

**Report** (under 400 words): commits; each item with its proof; each check with its exit code and test count; the
bounds under 50 ms you found; any expectation changed and why; what you left out; drift, gaps and followups.
