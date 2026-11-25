# Lane FX26: two fixes found overnight

**Model**: Sonnet (two closed items, each with its cause known). **Base**: `packages-seq` at the orchestrator's head
(`git reset --hard packages-seq` first). Branch `seq/fx26`. Beside you: FX25's worktree, stopped, holds uncommitted
work in `core/src/render/` ranges and lanes, `region-hover.ts`, `ext-rondo/` labels and the variant card: never edit
those files; if a fix needs one, stop and report.

**Why**: the loop's fix point after CL12: the open items of `fixes.md` that no lane holds and no decision waits on.

**Read first**: `.scaffolds/seq/execution/close-out/fixes.md`, the items "The legend's eye chunk lands after a test
file ends" and "A ring nobody follows narrows a hidden window"; `drift/ring-cmd-wheel-at-fit.md` (only its Evidence
and the shared part of its options: the window idle); `.scaffolds/seq/design/layouts.md` (Rings); the viewer spec's
Input, Legends and the views' zoom; the root `AGENTS.md`'s Tests section (files share their worker's registry).

**Items**, a commit at each:

1. **The eye chunk awaited where it starts** (`viewer/src/mount-legend.ts` `loadLegendEyes`, loaded on a legend's
   first mount): the class is a lazy chunk a test file starts and does not await, so whichever file mounts a legend
   first may end before the chunk lands (`EnvironmentTeardownError`, every test passing, exit 1). Fix the class: the
   tests that mount legends await the chunk (a shared helper the story tests already use, or the mount's own
   promise), not one file at a time. Proof: `pnpm test:unit packages/seq/` with `--sequence.shuffle.files` three
   times, exit 0, no teardown error; say how many files mount a legend.
2. **A ring nobody follows keeps its window** (`viewer/src/models/view.ts` `_zoomX` and the pan): zoom and pan idle on
   a ring section no navigator or follower uses, as `layouts.md` says; whether the wheel is then taken stays as it is
   now (the drift note decides it). Proof: tests that Cmd+wheel in and out and a plain wheel leave the section's
   viewport at the whole axis, and that a followed ring still moves.

**How**: host neutral; reuse before new helpers; module size; nothing per frame (allocation tiers green); each fix's
test fails on its parent (say so).

**Where**: a worktree on `seq/fx26`. Out of scope: FX25's files (above). Headers per the root `AGENTS.md`. Hot files
never staged: `package.json`s, `tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`, the work folder; `vitest.config.ts` and `vitest.setup.ts`
too (if item 1 needs them, stop and report). Never a bare `git stash`. No network.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf`; `pnpm -s -C packages/seq/core size` (every
row: headless 978.7 of 980, default 1,557.0 of 1,560, minimal 1,096.3 of 1,100); `npx tsc -b`; `pnpm lint`;
`pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report the
action and why; don't work around it).

**Report** (under 300 words): commits; each item with its proof and whether its test fails on the parent; the checks
with exit codes and counts; the size rows; drift, gaps and followups; the changelog text.
