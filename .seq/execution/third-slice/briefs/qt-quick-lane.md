# Lane QT: the quick lane's time (stages 3 and 4's gate)

**Goal**: the gate's time rule (`loop.md`, Gates; root `AGENTS.md`, Tests: three lanes): each touched folder's quick lane
under 2 s, the slowest files fixed when slow by waste, moved to the medium lane with a quick representative left when
slow by nature. Measured on `packages-seq` at the gate (`pnpm test:unit <folder>`, Duration): `packages/seq/core`
2.05 s, `packages/seq/storybook` 2.25 s (viewer 1.31, examples 1.58). The slowest files:

| ms | file |
|---|---|
| 2137 | `packages/seq/storybook/src/_stories/examples/_test/globins-stories.test.ts` |
| 2025 | `packages/seq/storybook/src/_stories/examples/_test/themes-stories.test.ts` |
| 1994 | `packages/seq/storybook/src/_stories/examples/_test/navigators-stories.test.ts` |
| 1868 | `packages/seq/core/src/theme/_test/ink.test.ts` |
| 1773 | `packages/seq/core/src/render/_test/rows-time.test.ts` (its quick half, the bounds are heavy) |
| 1335 | `packages/seq/storybook/src/_stories/examples/_test/navigators-mounted.test.ts` |
| 1237 | `packages/seq/core/src/render/_test/pairs-grids.test.ts` |

For each: find where the time goes (a fixture built per test or per `describe` instead of once per file; real timers
where fake ones would do; a full sweep where one representative would cover the quick lane), fix waste first, then move
what is slow by nature with `{ tags: ['medium'] }` (or `// @module-tag medium`) and keep a quick representative. Never
drop coverage: what moves runs in `pnpm test:medium`. Report each file's time before and after and what you did.

**Also measure, do not change the config**: the heavy lane's timing bounds fail only in the full parallel run and pass
alone (`rows-time`'s theme change, `pairs-grids-time`'s arcs hover; earlier `frame-time`). Measure whether running the
timing files without file parallelism (vitest's `fileParallelism: false` or `sequence` for the heavy mode, or a pool of
one for files that call `itWithin`) makes the full `pnpm test:heavy packages/seq/` pass three times in a row, and its
wall time. `vitest.config.ts` is a hot file: write the change you recommend in the report, the orchestrator applies it.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/qt-quick-lane
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/` (and each folder's Duration), `pnpm test:medium packages/seq/`, `npx tsc -b
packages/seq/core packages/seq/viewer packages/seq/examples`, `pnpm -C packages/seq/<core|viewer|examples>
typecheck:test`, `pnpm -C packages/seq/storybook typecheck`, `pnpm lint`, `pnpm format:check`. Before any heavy run,
`pgrep -fl "vitest.*heavy"`.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, every `package.json` and `tsconfig*.json`,
anything in `.seq/`.

**Stop when** each touched folder is under 2 s with the checks green, or a file cannot get there without dropping
coverage (report it), or a check fails twice after a fix.

**Report**: commits; each folder's Duration before and after; each file's; what moved to medium and its quick
representative; the heavy-lane measurement and the config change you recommend.
