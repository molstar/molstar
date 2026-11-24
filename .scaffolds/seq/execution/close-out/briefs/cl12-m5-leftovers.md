# Lane CL12: M5's leftovers

**Model**: strong (frame budgets: profile first, then the fix that holds across cases). **Base**: `packages-seq` at the
orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/cl12`. Beside you: none running; FX25's
worktree holds uncommitted work in `core/src/render/` ranges and lanes and the variant card: avoid those files; if you
must, stop and report.

**Why**: the fifth slice's measurements (M5) left three items over their targets or unanswered, carried to this plan
(`plan.md`, "M5's leftovers"; `packages/seq/TODO.md`, "M5's leftovers").

**Read first**: `.scaffolds/seq/execution/fifth-slice/measurements.md` (the method: the rAF queue, 20 warm and 300
measured frames, dpr 2, 1440 × 900; the rows "200 copies, split patch" and "by region (tabix bedGraph)", the gaps);
`.scaffolds/seq/execution/fourth-slice/measurements.md` (M4's method); `.scaffolds/seq/design/performance.md` (the
frame targets); `ux-foundations.md` (the wheel and the interaction matrix) and `ring-extension.md` (the ring's
gestures); the viewer spec's Input and History; the code the profiles name.

**Items**, a commit at each:

0. **Profile first**, each case in the M5 method, with the load said: the split patch on Examples › Complex › Views
   (Copies, then Split: click to drawn, 5 or more samples); Genomes › Views › Signal's zoom (300 frames; which frame is
   over and why: a tile level's first read, a layout, a GC); Cmd/Ctrl+wheel over the ring at the fit (what takes the
   event and what it does). For Cmd+wheel, state what the designs say it should do at the fit; if they don't say, stop
   that item with a drift note (`drift/`, the template in `.scaffolds/_method/drift-note.md`) with the options and a
   recommendation, and go on with the others.
1. **The split patch under 16 ms**, click to drawn (median and max): the work the patch does once, not per copy, or
   spread so the first frame stays in budget with the final picture following within the design's rule. Proof: a
   headless twin bounding its work (perf tier, allocation or operation count), and the browser numbers before and after.
2. **The Signal zoom's frames under 16 ms**: the over-budget frame's cause fixed by its class (a level's first build
   spread across frames, as FX1 did for rows, if that is the cause). Proof: a perf-tier case for the cause; the browser
   numbers before and after.
3. **Cmd/Ctrl+wheel over the ring at the fit**: does what item 0 found the designs say (or what the drift note's
   decision says, if the orchestrator relays one); never "taken but nothing". Proof: tests per direction.
4. **The records**: `performance.md`'s rows and the Storybook Performance page if a target's text moves; the viewer and
   ext-rondo specs for item 3. Proof: `designs-in-specs` green.

**How**: measure before changing code; fix the class, not the instance; host neutral (the viewer's frame scheduling,
never a global rAF); nothing per frame (allocation tiers green); module size (§module-size); browser numbers in the
Browser pane only when `uptime`'s load is under 3, the load said beside every number; a number read under load is
marked so, not used as a verdict.

**Where**: a worktree on `seq/cl12`. Out of scope: FX25's files (above). Headers per the root `AGENTS.md`. Hot files
never staged: `package.json`s, `tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`, the work folder (but a drift note goes in `drift/`, as the
loop allows). Never a bare `git stash`. No network.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` and `pnpm -s -C packages/seq/storybook
test:perf`; `pnpm -s -C packages/seq/core size` (every row; the minimal viewer has 0.3 kB left: new code lazy, or name
the row); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage
first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report the
action and why; don't work around it, not even by another tool).

**Report** (under 450 words): commits; item 0's profiles; each item with its proof and its numbers before and after
(with the load); the checks with exit codes and counts; the size rows; what to look at; drift, gaps and followups; the
changelog text.
