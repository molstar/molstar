# Lane QT3: the quick lane's time (after HC)

**Goal**: the gate's time rule again (`loop.md`, Gates; [testing](../../../design/testing.md): each folder's quick lane
under 2 s, a test file under 1 s, all of `packages/seq` under 30 s), as [QT2's brief](qt2-quick-lane.md) did. HC and NP
added to core, viewer and storybook. Measured on `packages-seq` after HC (`pnpm test:unit <folder>`, Duration, two runs
at load 3.1–4.8): `packages/seq/core` 2.40 / 2.69 s, `packages/seq/viewer` 2.47 / 2.22 s, `packages/seq/storybook`
2.63 / 2.08 s (examples 1.54 s). The slowest files (each alone in its folder's run):

| s | file |
|---|---|
| 0.86 | `packages/seq/core/src/render/_test/pairs-grids-time.test.ts` |
| 0.86 | `packages/seq/core/src/render/_test/navigator-time.test.ts` |
| 0.69 | `packages/seq/core/src/render/_test/pairs-grids.test.ts` |
| 0.60 | `packages/seq/core/src/render/_test/accessibility.test.ts` |
| 0.47 | `packages/seq/core/src/render/_test/rows-time.test.ts` |
| 0.46 | `packages/seq/core/src/render/_test/themes.test.ts` |
| 0.62 | `packages/seq/viewer/src/_test/changes.test.ts` |
| 0.38 | `packages/seq/viewer/src/_test/story.test.ts` |
| 0.37 | `packages/seq/viewer/src/_test/hover-card.test.ts` |
| 0.35 | `packages/seq/viewer/src/_test/drag-release.test.ts` |
| 0.80 | `packages/seq/storybook/src/_stories/features/text/_test/text-viewer.test.ts` |
| 0.80 | `packages/seq/storybook/src/_stories/features/themes/_test/themes-stories.test.ts` |
| 0.68 | `packages/seq/storybook/src/_stories/examples/_test/pairs-grids-stories.test.ts` |
| 0.60 | `packages/seq/storybook/src/_stories/features/story/_test/story-stories.test.ts` |

Re-measure first (each file alone and in its folder, with `uptime`; a load over 3 makes a number a guide only). For
each: fix waste first (a fixture or a viewer built per test or per `describe` instead of once per file; files read
again per test; real timers where fake ones do; a sweep over every case where one representative covers the quick
lane), then move what is slow by nature with `{ tags: ['medium'] }` or `// @module-tag medium`, keeping a quick
representative. The `*-time.test.ts` files: what in them is a wall-clock bound belongs under `itWithin` (heavy), its
unbounded half quick or medium. Never drop coverage: what moves runs in `pnpm test:medium`. Aim for room under the
guide (about 1.7 s a folder), so the next lane does not reopen it.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/qt3-quick-lane
<that commit>`. Nothing else runs in `packages/seq` beside you.

**Checks**, **headers**, **hot files** and **stop when**: as [QT's brief](qt-quick-lane.md) says (no heavy measurement
this time; `pnpm test:heavy packages/seq/` only if you move a bound into `itWithin`, after `pgrep -fl "vitest.*heavy"`).

**Report**: commits; each folder's Duration before and after (with the load); each file's; what moved to medium or
heavy and its quick representative; any waste fixed.
