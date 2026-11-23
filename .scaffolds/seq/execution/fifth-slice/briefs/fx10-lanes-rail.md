# Lane FX10: the lanes rail by gesture, and the keys

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-10`. Model:
Sonnet (a decided design). Follow `fx-fix-lane.md` (a commit per part, a test failing on its parent, a look in light
and dark).

**Why**: the author, 2026-11-22 (Features › Genome formats › GTF › TP53): "the scrollbar is nice, but it can't be
scrubbed by a gesture on the control itself, only by wheel". layouts.md's interaction matrix (line 650): "the wheel's
vertical delta over an overflowing lanes band scrolls its lanes; drag the slim rail; Page keys with the cursor in it".
FX9 built the rail (the rows rail's `placeUpright` geometry, `SlimRailPx` 12) and the wheel (`wheelLanes`, `framed` in
the rows behavior), not the gesture.

**Build**, reusing the rows rail's gesture (`rowsRailGesture`, `railUnitAt`, the `RowsRail` hit and `partCursor` in
`viewer/src/input.ts`), not a second copy:

1. **Drag the lanes rail's thumb**: the band's `bandTop` follows the pointer (grab, grabbing), clamped; Escape cancels
   to where it began (as the rows rail).
2. **Press on the rail outside the thumb**: pages toward the press (or jumps there, whichever the rows rail does: the
   same rule).
3. **Keys**: Page Up and Page Down with the cursor in an overflowing lanes band scroll it a page (the band's visible
   lanes less one); Home and End to the first and last lane if the rows rail has them.
4. **Hover**: the rail's thumb widens on hover as the rows rail's does; the cursor says so.

**Beside you**: LA (the bridge) writes in `viewer/src/input.ts` too (the bridge's pointer). Keep your change there a
dispatch to the rows behavior, small; put the work in the rows behavior and its tests. Stay out of `layout/bridge.ts`,
`pairs-ribbons.ts`, `models/view.ts`'s lock.

**Tests**: a drag on the rail scrolls the band and the stored `bandTop` (a headless twin on a fake clock); a press
outside the thumb pages; Page keys; Escape cancels; the drag allocates nothing warm (both perf tiers: a scenario "gene
models, rail dragged"). Each fails on its parent.

**Where**: a worktree on `seq/fx-10`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`.scaffolds/seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s, `tsconfig*.json`): list edits. Stop any
Storybook server you start and close any browser tab you open.

**Checks** (with numbers): the quick lane (under 2 s a touched folder), `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, both perf tiers, `pnpm -s -C packages/seq/core size` (headless 923.3 of 930, default 1,464.2 of
1,470, minimal 1,039.3 of 1,040 kB: a bound passed stops the lane with a report), the `.scaffolds/seq/_checks`
(features-planned with `--plan .scaffolds/seq/execution/fifth-slice/plan.md`).

**Report**: under 250 words: commits, the tests that failed on the parent, sizes, perf, the look, the
`.v6/changelogs/` text.
