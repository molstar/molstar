# Lane FX31: the reversed pan's leak control fails at GATE 2

**Model**: Sonnet (one perf case, bisected, its cause fixed). **Base**: `packages-seq` at the orchestrator's head
(`git reset --hard packages-seq` first). Branch `seq/fx31`. No lane beside you.

**Why**: GATE 2 (the orchestrator, 2026-11-24, after FX29, FX30 and FX28 merged). `pnpm -s -C packages/seq/core
test:perf` failed one case, and the file alone failed it twice more at load 2.2–2.6 (so not a flake: a fix lane's):

```
FAIL render/_test/frame.quick.ts > … > reversed, 150,000 pieces: pan (layout and draw)
AssertionError: reversed, 150,000 pieces pan +leak: windows 391345.32 / 391933.52 / 391797.00 B (200 steps each),
warmed 200 steps (capped), 240 GC(s) accounted, code 9.37 B: expected 2006.2 to be greater than 3072
```

The case's own reading (bytes a frame, under its `LayoutBound`) passes; the planted-leak control
(`test-harness/src/quick.ts`, `quickCase`: the measured delta must exceed 0.75 of `controlBytes`) reads 2.0–2.4 kB of
a 4 kB leak. It passed at GATE 1 (`b6d652e27`) and in FX28's lane (before FX29 and FX30 merged).

**Read first**: `test-harness/src/quick.ts` (`quickCase`, `scaledBytes`, the control's rule and its tests),
`packages/seq/core/src/render/_test/frame.quick.ts` (the scenarios, `LayoutBound`), the merges since GATE 1 on
`packages-seq` (`git log --first-parent b6d652e27..packages-seq`).

**Items**, a commit at each:

0. **Bisect**: the first-parent merges since GATE 1, this case alone (`npx vitest run -c vitest.perf.config.ts
   render/_test/frame.quick.ts -t 'reversed, 150,000 pieces: pan'` in `packages/seq/core`), at load under 3 (`uptime`
   first; wait if higher). Say the merge, then the commit inside it. Proof: the readings at the last good and the first
   bad commit.
1. **The cause, fixed**: if the layout's allocation grew or became noisier (more collections a window), that is the
   regression: fix it in the code that grew, never by loosening a bound or the control's rule. If the control itself
   mis-scales for a case of ~390 kB a step (its `scaledBytes`), fix the harness's rule with its test from given
   readings, and say why the case's reading is unchanged. Proof: the case green three times alone; the full core
   `test:perf` green.
2. **The records**: `testing.md` or the harness's comment if the control's rule moves. Proof: the scaffold checks.

**How**: allocation bounds are hard (never raised to pass); nothing per warm frame; host neutral.

**Where**: a worktree on `seq/fx31`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network. Do not run the heavy lane.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/
test-harness/`; `pnpm -s -C packages/seq/core test:perf` and `pnpm -s -C packages/seq/storybook test:perf` (load
under 3); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage
first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 300 words): the bisect's result and the cause; commits; each item with its proof; the checks with
exit codes, counts and the load; followups; the changelog text, if any.
