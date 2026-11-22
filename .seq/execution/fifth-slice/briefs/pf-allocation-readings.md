# Lane PF: allocation readings that depend on the process

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/pf`. Model:
strong (measurement, V8, the harness). HC is in hover paths across layers, LA3 in the grid layer: change their code
only where a real allocation is found, and list each touch.

**Why**: `fixes.md` Open, "Allocation readings that depend on the process, not the code" (read it whole). The
allocation tiers are a merge gate (the author: "the important parts are the performance and low allocations"); three
cases flip by what ran before them in their file, so the gate blocks merges at random or hides a real allocation:

- storybook `_stories/features/layouts/_test/section-scroll.quick.ts`, the wheel step: 0 B alone, 40 B and 80 B in a
  full-file run (FX14's two-case probe, wheel over an empty lane and over a gene, is at
  `fx14-two-case-probe.ts` in the orchestrator's scratchpad (ask for it if you need it));
- storybook `chromosome` "between two levels" pan: 55,9xx–56,1xx B against its 56,000 bound;
- core `render/_test/navigator-blocks.quick.ts`, the 3,000-contig leak control: ~112 KB, failing in alternate runs.

**Task**:

1. For each case, find the cause with evidence (a heap profile or sampling of allocations by site, `--trace-gc`, the
   harness's own counters): a real allocation in our code (then fix it, a test at 0 B that fails on its parent), V8
   state from earlier cases (deopts, inline caches, lazy feedback vectors), or the harness's warm-up, window or GC
   accounting.
2. If the harness is the cause: the harness is `test-harness/` (`@molstar/test-harness`, shared by every package). Do
   not change it in this lane: write the change you would make (the file, the diff in prose, why, what it changes for
   the other packages' tiers) in the report, and make the seq cases robust on their side if that is honest (each case
   warmed in its own order, or isolated by file) without loosening a bound.
3. The chromosome case: if it is real growth, say from which commit (bisect over `packages-seq`'s merges since the
   bound was set) and what allocates; if noise, show its distribution over 10 runs and what makes it vary.

No bound is raised or loosened. Timing bounds are untouched.

**Where**: a worktree on `seq/pf`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`). Anything in `@molstar/core`, `@molstar/io` or `test-harness/` stops at a report. The machine runs
other lanes: measure each case in a fresh process, several times, and say the load average beside each reading.

**Checks**: the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`, each
`typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers to completion three times
each (every reading listed), `pnpm test:medium packages/seq/`.

**Report**: under 400 words: per case, the cause with its evidence, the fix or the proposed harness change; commits;
the three full-tier runs' readings of the three cases; anything else found.
