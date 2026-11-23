# Lane FX20: the shared harness's leak control scaled to the case

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx20`. Model:
Sonnet. **Scope**: `test-harness/` (`@molstar/test-harness`, shared by every package; the author approved this change,
2026-11-23, stage 3's gate) and seq's tests that use its own `scaledControl`.

**Why** (PF, fifth slice): `quickCase` measures the case, then `leaking(case)`, as two medians; the 88-byte control is
below the harness's resolution for cases allocating KB a step (windows wander 20–250 B a step with hundreds of GCs), so
those cases carry `control: false` and seq's own `packages/seq/core/src/render/_test/scaled-control.ts` (a 4 KB array a
step, checked ≥ 3 KB).

**Task**: in `test-harness/src/quick.ts` (read it and its tests first): when a case reads above ~1 KB a step or its
windows spread above ~8 B, size the control's object to ~1% of the case (at least 4 KB) and check 75% of it; or, if
cleaner and as cheap, alternate case and control windows on the same warm state and take the median of the pairs'
differences. Cases that allocate nothing keep today's 88 B control unchanged. Then seq's `scaledControl` and its
`control: false` cases (navigator-blocks, bridge relayout, categorical grid pan, 2OCJ hover, and the storybook cases
PF listed: a bridge ribbon hover, a pan through a lock, a connector, linked pairs) use the harness's control; remove
`scaled-control.ts`. The harness's own tests: a KB-a-step case with a real leak scaled to it is caught, one without
passes, a zero-allocation case unchanged.

**Where**: a worktree on `seq/fx20`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`). Nothing in `@molstar/core` or `@molstar/io`. No bound
loosened.

**Checks**: `pnpm test:unit test-harness/ packages/seq/ cli/seqx-render/`, `npx tsc -b`, `pnpm lint`,
`pnpm format:check`, both seq perf tiers twice to completion (list each moved case's reading), and every other package's
opt-in perf tier that uses `quickCase` once (find them; list each result) so no other package breaks.

**Report**: under 200 words: the design chosen and why, commits, the readings, other packages' results, the changelog
line (`.v6/changelogs/`, for the orchestrator).
