# Lane FX3: the positions rail's place, and the open fixes

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-3`. Follow
`fx-fix-lane.md` (one commit per item, a test failing on its parent, a look in light and dark). Read core's and the
viewer's `_spec` (Render: Navigators, Rings, Blocks; viewer Configuration) and `.seq/design/navigators.md` ("The
positions rail beside a declared navigator, and its place", and the Configuration row) first.

**Items** (`fixes.md`, Open), in this order:

1. **The positions rail beside a declared navigator, and its place** (designed 2026-11-21, navigators.md):
   `navigator.positionsRail` (`ViewerConfig`: `always` default, `unlessNavigator`, `none`) and `navigator.railPlace`
   (`ruler` default, `bottom`, `navigator`). The Globins rows family story (`seq-next-examples-globins-rows--family`)
   then shows a horizontal rail beside its vertical overview; add a story control or a second story showing each place.
   Twin tests: the rail present with a declared overview under `always`, absent under `unlessNavigator`; its row in
   each place; a drag on it pans and the overview's box follows in the same frame. Viewer and core `_spec` updated
   (the design's text is the rule).
2. **A ring's own hover draws no mark** (three-thousand-contigs: hovering the ring shows the card but no position box on
   the ring; hovering the row draws it there). The hovered view draws its own hover mark. Twin test on a ring and a row.
3. **An assembly `sequence` override naming no sequence is dropped silently** (GC1's find): raise
   `seq.resolve.assembly.*` as for an alias naming nothing; a resolve test.
4. **Block rulers unnumbered at the fit** (the-hosts-complex: chain blocks ~180 px show ticks without numbers in a real
   browser; headless draws "100"/"200"): look in the browser first; if it reproduces, find the text-metrics path that
   differs and fix it with a test that fakes the browser's metrics; if it does not reproduce, say so with a screenshot.

**Where**: a worktree on `seq/fx-3`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files
never staged (`.seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s): list the edits. Public API changes
allowed in `packages/seq`; anything else stops the item.

**Checks**: `pnpm test:unit packages/seq/`, `npx tsc -b`, `pnpm -s -C packages/seq/{core,viewer,examples}
typecheck:test`, `pnpm -C packages/seq/storybook typecheck`, `pnpm lint`, `pnpm format:check`, both perf tiers,
`pnpm -s -C packages/seq/core size`, `pnpm test:heavy packages/seq/`, the four `.seq/_checks`. Run them all before the
report and give numbers.

**Report**: under 250 words: commits per item; checks with numbers; the stories to look at; anything left and why; the
`.v6/changelogs/` text.
