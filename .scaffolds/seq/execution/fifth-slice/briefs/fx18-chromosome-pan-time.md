# Lane FX18: the chromosome pan between two levels, back under its time bound

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx18`. Model:
Sonnet. No other lane is running (a gate is open): measure on a quiet machine (`uptime` under 3; say the load beside
every reading).

**Why** (stage 2's and 3's gate, 2026-11-23): `pnpm test:heavy packages/seq/examples/src/_test/chromosome.test.ts`
fails "a pan frame within a bound of the p53 gene view's" at "between two levels" (ratio 1.9, `chromosome.test.ts:498`):
best ~21.7–22.1 ms against 1.9× ~11.0–11.6 ms, in every run, alone, at load ~2.2. Timing bounds are gates: never
loosen it.

**Task**:

1. Find the merge that pushed it over: bisect the first-parent merges on `packages-seq` since the stage-1 gate (the
   suspects, in order: LA4 `seq(fifth-slice): merge LA4` (tracks over an assembly sized by every block's lanes), FX14
   (close-sharing clusters), HC (`drawOwnHover`, the sweep's fixes), VL/VP (the values look), FX16, FX17, TT (card
   providers, tooltips' door)); a temporary worktree per commit, `pnpm install --offline`, the one test, three runs
   each. Report the ratio per merge.
2. Profile the guilty frame (`--cpu-prof` or the harness's timing) and fix the cause in the code (work per frame that
   should be cached at layout or per settle; a per-block loop that should run once), keeping its behaviour and its
   tests. A test that would catch it again (a count of the work per pan frame, quick lane) failing on its parent.
3. Rerun the heavy file three times alone: every case under its ratio.

**Where**: a worktree on `seq/fx18`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`). Anything in `@molstar/core` or `@molstar/io` stops the lane.
Remove the temporary worktrees you create.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, `pnpm lint`, `pnpm format:check`, the core perf tier, the heavy chromosome file three times.

**Report**: under 200 words: the guilty merge with the per-merge ratios, the cause and the fix, the three heavy runs'
readings with the load, commits.
