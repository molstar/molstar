# Lane QT2: the quick lane's time (stage 5's gate)

**Goal**: the gate's time rule again (`loop.md`, Gates; [testing](../../../design/testing.md): each folder's quick lane
under 2 s, a test file under 1 s, all of `packages/seq` under 30 s), as [QT's brief](qt-quick-lane.md) did after stages
3 and 4. Measured on `packages-seq` at stage 5's gate (`pnpm test:unit <folder>`, Duration): `packages/seq/examples`
2.04 s, `packages/seq/storybook` 2.45 s (core 1.91, viewer 1.86; 6.37 s in all). The files over 1 s (their summed
time in a full run):

| s | file |
|---|---|
| 2.30 | `packages/seq/storybook/src/_stories/features/themes/_test/themes-stories.test.ts` |
| 1.54 | `packages/seq/examples/src/_test/navigators.test.ts` |
| 1.45 | `packages/seq/examples/src/_test/views.test.ts` |
| 1.02 | `packages/seq/examples/src/_test/plasmid.test.ts` |
| 1.00 | `packages/seq/storybook/src/_stories/examples/_test/complex-stories.test.ts` |

Re-measure first (each file alone and in its folder); the datasets lanes (DG, DC, DA, CY's views over the files) added
most of this. For each: fix waste first (a fixture or a viewer built per test or per `describe` instead of once per
file; files read again per test; real timers where fake ones do; a sweep over every view where one representative
covers the quick lane), then move what is slow by nature with `{ tags: ['medium'] }` or `// @module-tag medium`,
keeping a quick representative. Never drop coverage: what moves runs in `pnpm test:medium`. Core and viewer are near
the guide: keep them under it.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/qt2-quick-lane
<that commit>`. Nothing else runs in `packages/seq` beside you but the novice pass (read-only, in a browser).

**Checks**, **headers**, **hot files** and **stop when**: as QT's brief says (no heavy measurement this time).

**Report**: commits; each folder's Duration before and after; each file's; what moved to medium and its quick
representative.
