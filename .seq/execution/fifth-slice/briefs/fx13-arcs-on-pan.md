# Lane FX13: arcs that hold their shape while the view moves

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-13`. Model:
strong (arc geometry and the pairs layer's culling and levels; the previous Seq* had the same bugs and fixed them).

**Why**: the author, 2026-11-22 (Examples › Complex, the Contacts (≤ 4 Å) track over CDK2 (A) and cyclin A (B), the
rail's window dragged): "the arc rendering doesn't really work when moving the rail. similar issues were present in
packages/sequences too. bad these came back again." The screenshots: at one window every arc gone, only grey ticks at
the band's foot where the partners were; at the next, arcs back; an arc whose partner is past the right edge drawn as
a long line rising off the band; the set of arcs changing as the window moves a little.

**The rules, as `packages/sequences` settled them** (read `packages/sequences/viewer/src/draw/arcs.ts`,
`renderer-helpers.ts` (`Arcs.frame`, `ARC_MARGIN` 0.25, `ARC_PAST` 4) and its tests `viewer/_test/arcs-off-window.test.ts`,
`arcs-continuous.test.ts`, `arcs-pan.test.ts`; port the rules, not the code):

1. **A partner past the edge**: the arc is drawn over its virtual end, the partner taken at most `ARC_MARGIN` of the
   window past the edge; never flat, never out of the band; a pair spanning the whole window is the widest arc over
   it.
2. **Continuous under a pan**: panned a px at a time, across the moments a partner enters or leaves the window and the
   widest arc enters or leaves the view, every arc's points move a px or so, never jump; the arcs' scale (the widest
   span in view) does not snap as the widest arc comes and goes.
3. **A pan's frame equals a full frame**: the scale is the whole view's, never a strip's.
4. **Items, not ticks, while items fit**: a pair whose one end is in view and whose other end is off it stays an arc
   (rule 1), never a tick at the band's foot; the switch to bins (past the gate) is by the pairs in the reachable
   region at the level, not by how many have an end in the window, and does not flicker as the window moves.

Find which of these the new pairs layer (`core/src/render/pairs*.ts`; not `pairs-ribbons*`, LA2's) breaks, and fix it.

**Tests** (port the three test files' cases to the new layer's fixtures: the complex's contacts, a junction track):
a partner off-window (virtual ends, no flat arc, inside the band), continuity (sample each arc's points across a pan
a px at a time: no jump over 2 px), a pan's frame against a fresh full frame, items never turned to ticks while their
count fits as the window moves; canvas and SVG; each failing on its parent. Both perf tiers (a warm pan 0 B new).

**Beside you**: LA2 (bridges: `pairs-ribbons*`, `pairs-region.ts`, resolve, the store), FX12 (empty wells:
`tracks*.ts`, `blocks.ts`, `pairs.ts:~425` the ring's past-cap chords), VL (the values layer). Keep out of their
files; if `pairs.ts` must change, keep it to the arc path and list the hunk.

**Where**: a worktree on `seq/fx-13`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged: list them. Stop any Storybook server you start and close any browser tab you open. Sizes (headless 940,
default 1,480, minimal 1,060 kB): read the size rows for FAIL (the exit code is 0 on a failing row).

**Checks**: the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`, each
`typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers to completion, size,
`pnpm test:medium packages/seq/`, the `.seq/_checks`.

**Report**: under 300 words: which rules were broken and why, commits, the tests that failed on the parent, the look
(Examples › Complex contacts, the rail dragged across), sizes.
