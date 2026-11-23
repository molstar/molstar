# Lane LA: layers across sections (the bridge first)

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/la`. Model:
strong (a new layout scope, the hit, the viewer's lock). **FX9** (track sizing, the tree's way back) is finished but not
merged: its files are `render/layout/tracks.ts`, the rows behavior, the ranges layer and `known.ts`, `narrow.ts`,
`header-extras.ts`, `viewer.ts`'s `renderSvg`. Touch them only where the bridge must (a stacked layout's entry, export);
list each such touch in the report, as the orchestrator merges FX9 into your branch when the author decides its sizes.

**Task**: [across-sections](../../../design/across-sections.md)'s Consequences, in this order, a green commit at each
step; **stop after step 4** if the lane runs long and report what is left (steps 5–6 become LA2):

1. **Resolve**: the `bridge` kind as GC1 declared it (`a`, `b`, `label`, `height`, `lock`, `hidden`; children
   `layer`); a `pairs` layer's `data` a segments `map`; `pairs.shape`; `ResolvedView.bridges` and
   `ResolvedView.stack`; the issues the design names (`seq.lang.bridge.sides`, `seq.view.bridge.*`,
   `seq.render.bridge.unplaced`) in the catalog. **No language bump**: the params are declared; a missing declaration
   stops the lane with a report.
2. **Layout and render**: `layout/bridge.ts` (placement between neighbouring sections, its band, its ends);
   `LayerType.scope: 'bridge'` and `BridgeContext`; `pairs-ribbons.ts` (ribbons and stubs, the twist by GA's rule,
   colour by the a end's block, the per-side reads bounded by `maxTilesPerLayer`, the drawn-once rule over items and
   bins, emphasis); the closed-form hit and `Hit.bridge`; export (a bridge a group with `data-bridge`, `data-a`,
   `data-b`; ribbons as cubic paths). Lazy: nothing of it in a view without a bridge (a size row for the chunk,
   about 6 kB).
3. **Viewer**: the bridge's pointer (drag and wheel), `lock` in the `ViewModel` beside follows and groups,
   `bridge.align` (`navigate`), `bridge.selectAcross` (`select`), the lock toggle, the cursor across a bridge, the
   cards, ARIA; the bindings checked (`Shift+A`).
4. **Examples and stories**: peach and grape with a bridge (an inversion, a translocation, a gap), one genome against
   itself with `mirror`, three genomes chained (generated in `seqx-examples`, small, deterministic on `PCG`);
   `Features/Across sections/*` one story per Default, light and dark, canvas and SVG, each with its headless twin on
   a fake clock (tiles arriving per side); a Design page "Across sections" as the other Design pages.
5. **Grids over two sections**: a grid section's `follow_y` a second section of blocks (a dot plot of two genomes, a
   Hi-C between two chromosome sections); categorical grid axes (`of`/`of_y` a rows system, `member_order`,
   `member_order_y`; a kinship grid beside the genotype matrix).
6. **The liftover** last (EC-212): generated, its truths computed in code.

**Tests** (the design's Tests line): the closed-form hit equals the drawn polygon's topmost cover; the twist by
geometry equals GA's three-direction rule over every combination of block and section direction; each pair drawn
exactly once across two sides at two levels (windows overlapping or not, a tie in scale); a bridge's reads equal the
whole read clipped to either window; `lock` round trips; `align` lands on the hull it names; every ribbon's ends at its
data's positions (at the fit, across a tile edge, a block gap, reversed blocks); EC-129, EC-169, EC-206, EC-209,
EC-212, EC-218, EC-223, EC-235 in titles where they apply. **Performance**: a bridge frame of 4,096 ribbons under 8 ms
(heavy `itWithin`; 1,024 quick), its hover under 2 ms at 0 bytes warm (an allocation case), a pan of one side relaying
that side and the bridge only (counted).

**Where**: a worktree on `seq/la`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each
green step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`): list them. The `_spec` sections the design lists
change in the same commits; other designs' amendments listed for the orchestrator (they are in `.seq/`). Anything in
`@molstar/core` or `@molstar/io` stops the lane. A file past 800 lines is a question, past 1500 a failure. **Sizes**
(reading 260, headless 920 — at 919.9 —, default 1,460 — at 1,459.6 —, minimal 1,040 kB; the author is deciding a
raise): the bridge is lazy; any growth of an eager bundle is reported with its cause; do not edit the bounds.

**Checks**: the quick lane (under 2 s a touched folder), `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf
tiers, `pnpm -s -C packages/seq/core size`, `pnpm test:medium packages/seq/`, `pnpm test:heavy packages/seq/` when load
is under 6, the `.seq/_checks` (features-planned with `--plan .seq/execution/fifth-slice/plan.md`). A fix's test fails
on its parent.

**Report**: under 400 words: commits; the steps done and left; checks with numbers; story ids; sizes (the chunk's row);
the public surface; the FX9 files touched; designs to amend; hot-file needs; the `.v6/changelogs/` text.
